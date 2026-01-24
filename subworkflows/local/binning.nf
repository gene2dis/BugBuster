/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BINNING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenomic binning using multiple tools and refinement with MetaWRAP
    
    DESCRIPTION:
        This subworkflow performs mode-agnostic metagenomic binning. It accepts
        contigs and BAM files (regardless of assembly or co-assembly origin) and
        runs a unified binning pipeline.
    
    INPUTS:
        contigs_and_bam: [meta, contigs, bam] - Contigs with aligned reads (BAM)
        checkm2_db: CheckM2 database for quality assessment
        gtdbtk_db: GTDB-Tk database for taxonomic classification
        reads: [meta, reads] - Original reads (optional, only for co-assembly bin coverage)
    
    OUTPUTS:
        refined_bins: [meta, bins] - MetaWRAP-refined bins
        versions: Tool version information
    
    WORKFLOW:
        1. Calculate contig depth from BAM files
        2. Run three binning tools in parallel (MetaBAT2, SemiBin, COMEBin)
        3. Refine bins with MetaWRAP
        4. Assess quality with CheckM2
        5. Classify taxonomy with GTDB-Tk
        6. Generate reports (mode-specific: simple reports for assembly, 
           coverage analysis for co-assembly)
----------------------------------------------------------------------------------------
*/

// Modules - unified for both per-sample and co-assembly modes
include { CALCULATE_DEPTH            } from '../../modules/local/calculate_depth/main'
include { METABAT2                   } from '../../modules/local/metabat2/main'
include { SEMIBIN                    } from '../../modules/local/semibin/main'
include { COMEBIN                    } from '../../modules/local/comebin/main'
include { METAWRAP                   } from '../../modules/local/metawrap/main'
include { CHECKM2                    } from '../../modules/local/checkm2/main'
include { GTDB_TK                    } from '../../modules/local/gtdb-tk/main'
include { BOWTIE2_SAMTOOLS_DEPTH     } from '../../modules/local/bowtie2_samtools/main'
include { BEDTOOLS                   } from '../../modules/local/bedtools/main'
include { BIN_QUALITY_REPORT         } from '../../modules/local/bin_quality_report/main'
include { BIN_TAX_REPORT             } from '../../modules/local/bin_tax_report/main'
include { BIN_SUMMARY                } from '../../modules/local/bin_summary/main'

workflow BINNING {
    take:
    contigs_and_bam // channel: [ val(meta), path(contigs), path(bam) ]
    checkm2_db      // channel: path(checkm2_db)
    gtdbtk_db       // channel: path(gtdbtk_db)
    reads           // channel: [ val(meta), [ reads ] ] - optional, only for co-assembly bin coverage

    main:
    ch_versions = channel.empty()

    //
    // Unified binning workflow - mode-agnostic
    //
    
    // Calculate depth from BAM files
    ch_depth = CALCULATE_DEPTH(contigs_and_bam)
    ch_versions = ch_versions.mix(CALCULATE_DEPTH.out.versions.first())
    
    // Run binning tools in parallel
    ch_metabat2 = METABAT2(ch_depth.depth)
    ch_semibin = SEMIBIN(contigs_and_bam)
    ch_comebin = COMEBIN(contigs_and_bam)
    
    // Collect versions
    ch_versions = ch_versions.mix(METABAT2.out.versions.first())
    ch_versions = ch_versions.mix(SEMIBIN.out.versions.first())
    ch_versions = ch_versions.mix(COMEBIN.out.versions.first())

    // Combine bins from all binners for refinement
    ch_all_bins = ch_metabat2.bins
        .join(ch_semibin.bins)
        .join(ch_comebin.bins)
    
    // Refine bins with MetaWRAP
    ch_metawrap = METAWRAP(ch_all_bins)
    ch_refined_bins = ch_metawrap.bins
    ch_versions = ch_versions.mix(METAWRAP.out.versions.first())

    // Quality assessment with CheckM2
    ch_checkm = CHECKM2(
        ch_all_bins.join(ch_metawrap.bins).combine(checkm2_db)
    )
    ch_versions = ch_versions.mix(CHECKM2.out.versions.first())

    // Taxonomic classification with GTDB-TK
    ch_gtdb_tk = GTDB_TK(ch_metawrap.bins.combine(gtdbtk_db))
    ch_versions = ch_versions.mix(GTDB_TK.out.versions.first())

    //
    // Generate reports - mode-specific handling
    //
    if ( params.assembly_mode == "assembly" ) {
        // Per-sample mode: simple quality and taxonomy reports
        BIN_QUALITY_REPORT(ch_checkm.all_reports.map { _meta, reports -> reports }.collect())
        BIN_TAX_REPORT(ch_gtdb_tk.report.map { _meta, reports -> reports }.collect())
    } else if ( params.assembly_mode == "coassembly" ) {
        // Co-assembly mode: additional bin coverage analysis and summary report
        ch_bin_depth = BOWTIE2_SAMTOOLS_DEPTH(
            reads.combine(ch_metawrap.bins.map { _meta, bins -> bins })
        )
        ch_bin_cov = BEDTOOLS(ch_bin_depth)

        BIN_SUMMARY(
            ch_bin_cov.collect()
                .combine(ch_gtdb_tk.report.map { _meta, reports -> reports }.collect())
                .combine(ch_checkm.metawrap_report.map { _meta, reports -> reports }.collect())
        )
    }

    emit:
    refined_bins = ch_refined_bins    // channel: [ val(meta), path(bins) ] or path(bins)
    versions     = ch_versions        // channel: path(versions.yml)
}
