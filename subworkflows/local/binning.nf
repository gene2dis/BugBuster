/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BINNING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenomic binning using multiple tools and refinement with MetaWRAP
----------------------------------------------------------------------------------------
*/

// Unified modules - handle both per-sample and co-assembly modes
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
    reads           // channel: [ val(meta), [ reads ] ]
    contigs         // channel: [ val(meta), path(reads), path(contigs) ] or path(contigs)
    bam             // channel: path(bam) or [ val(meta), path(bam) ]
    checkm2_db      // channel: path(checkm2_db)
    gtdbtk_db       // channel: path(gtdbtk_db)

    main:
    ch_versions = Channel.empty()
    ch_refined_bins = Channel.empty()

    //
    // Per-sample binning mode
    //
    if ( params.assembly_mode == "assembly" ) {
        // Calculate depth
        ch_depth = CALCULATE_DEPTH(bam)
        
        // Run binning tools
        ch_metabat2 = METABAT2(ch_depth.depth)
        ch_semibin = SEMIBIN(bam)
        ch_comebin = COMEBIN(bam)
        
        // Collect versions
        ch_versions = ch_versions.mix(METABAT2.out.versions.first())
        ch_versions = ch_versions.mix(SEMIBIN.out.versions.first())
        ch_versions = ch_versions.mix(COMEBIN.out.versions.first())

        // Combine bins for refinement
        ch_all_bins = ch_metabat2.bins.join(ch_semibin.bins).join(ch_comebin.bins)
        
        // Refine bins with MetaWRAP
        ch_metawrap = METAWRAP(ch_all_bins)
        ch_refined_bins = ch_metawrap.bins
        
        // Collect versions
        ch_versions = ch_versions.mix(METAWRAP.out.versions.first())

        // Quality assessment with CheckM2
        ch_checkm = CHECKM2(
            ch_all_bins.join(ch_metawrap.bins).combine(checkm2_db)
        )
        
        // Collect versions
        ch_versions = ch_versions.mix(CHECKM2.out.versions.first())

        // Taxonomic classification with GTDB-TK
        ch_gtdb_tk = GTDB_TK(ch_metawrap.bins.combine(gtdbtk_db))
        
        // Collect versions
        ch_versions = ch_versions.mix(GTDB_TK.out.versions.first())

        // Generate reports
        BIN_QUALITY_REPORT(ch_checkm.all_reports.map { meta, reports -> reports }.collect())
        BIN_TAX_REPORT(ch_gtdb_tk.report.map { meta, reports -> reports }.collect())
    }

    //
    // Co-assembly binning mode
    //
    if ( params.assembly_mode == "coassembly" ) {
        // Create coassembly meta
        def coassembly_meta = [id: 'coassembly']
        
        // Prepare inputs with coassembly meta
        ch_bam_collected = bam.collect()
        ch_contigs_with_meta = channel.of(coassembly_meta)
            .combine(contigs)
            .combine(ch_bam_collected)
        
        // Calculate depth
        ch_depth_co = CALCULATE_DEPTH(ch_contigs_with_meta)
        
        // Run binning tools with coassembly meta
        ch_metabat2_co = METABAT2(ch_depth_co.depth)
        ch_semibin_co = SEMIBIN(ch_contigs_with_meta)
        ch_comebin_co = COMEBIN(ch_contigs_with_meta)
        
        // Collect versions
        ch_versions = ch_versions.mix(METABAT2.out.versions.first())
        ch_versions = ch_versions.mix(SEMIBIN.out.versions.first())
        ch_versions = ch_versions.mix(COMEBIN.out.versions.first())

        // Combine bins for refinement
        ch_all_bins_co = ch_metabat2_co.bins.join(ch_semibin_co.bins).join(ch_comebin_co.bins)
        
        // Refine bins with MetaWRAP
        ch_metawrap_co = METAWRAP(ch_all_bins_co)
        ch_refined_bins = ch_metawrap_co.bins
        
        // Collect versions
        ch_versions = ch_versions.mix(METAWRAP.out.versions.first())

        // Quality assessment
        ch_checkm_co = CHECKM2(
            ch_all_bins_co.join(ch_metawrap_co.bins).combine(checkm2_db)
        )
        
        // Collect versions
        ch_versions = ch_versions.mix(CHECKM2.out.versions.first())

        // Taxonomic classification
        ch_gtdb_tk_co = GTDB_TK(ch_metawrap_co.bins.combine(gtdbtk_db))
        
        // Collect versions
        ch_versions = ch_versions.mix(GTDB_TK.out.versions.first())

        // Calculate bin coverage
        ch_bin_depth = BOWTIE2_SAMTOOLS_DEPTH(reads.combine(ch_metawrap_co.bins.map { meta, bins -> bins }))
        ch_bin_cov = BEDTOOLS(ch_bin_depth)

        // Generate summary report
        BIN_SUMMARY(
            ch_bin_cov.collect()
                .combine(ch_gtdb_tk_co.report.map { meta, reports -> reports }.collect())
                .combine(ch_checkm_co.metawrap_report.map { meta, reports -> reports }.collect())
        )
    }

    emit:
    refined_bins = ch_refined_bins    // channel: [ val(meta), path(bins) ] or path(bins)
    versions     = ch_versions        // channel: path(versions.yml)
}
