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
include { CHECKM2_BATCH              } from '../../modules/local/checkm2_batch/main'
include { GTDB_TK_BATCH              } from '../../modules/local/gtdb_tk_batch/main'
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

    // Collect all bins for batched processing
    // Combine bins from all binners for each sample
    ch_all_sample_bins = ch_metabat2.bins
        .join(ch_semibin.bins)
        .join(ch_comebin.bins)
        .join(ch_metawrap.bins)
        .map { meta, metabat, semibin, comebin, metawrap ->
            // Stage bins with sample-specific directories
            [
                meta,
                metabat, semibin, comebin, metawrap
            ]
        }
        .collect()
        .map { items ->
            // Extract metadata list and organize bins by sample and binner
            def meta_list = items.collect { item -> item[0] }
            def all_bins = []
            
            items.each { item ->
                def meta = item[0]
                def metabat = item[1]
                def semibin = item[2]
                def comebin = item[3]
                def metawrap = item[4]
                
                // Add bins with sample-specific staging paths
                all_bins.add(["${meta.id}_metabat", metabat])
                all_bins.add(["${meta.id}_semibin", semibin])
                all_bins.add(["${meta.id}_comebin", comebin])
                all_bins.add(["${meta.id}_metawrap", metawrap])
            }
            
            [meta_list, all_bins.collect { bin_tuple -> bin_tuple[1] }]
        }
    
    // Batched quality assessment with CheckM2
    ch_checkm_batch = CHECKM2_BATCH(ch_all_sample_bins.combine(checkm2_db))
    ch_versions = ch_versions.mix(CHECKM2_BATCH.out.versions.first())
    
    // Transform batched output back to per-sample format
    ch_checkm_all_reports = ch_checkm_batch.all_reports
        .flatMap { meta_list, reports ->
            // Group reports by sample
            def sample_reports = [:]
            reports.each { report ->
                def sample_id = report.name.split('_')[0]
                if (!sample_reports.containsKey(sample_id)) {
                    sample_reports[sample_id] = []
                }
                sample_reports[sample_id].add(report)
            }
            // Emit per-sample tuples
            sample_reports.collect { sample_id, report_list ->
                def meta = meta_list.find { m -> m.id == sample_id }
                [meta, report_list]
            }
        }
    
    ch_checkm_metawrap = ch_checkm_batch.metawrap_report
        .flatMap { meta_list, reports ->
            // Group reports by sample
            def sample_reports = [:]
            reports.each { report ->
                def sample_id = report.name.split('_')[0]
                if (!sample_reports.containsKey(sample_id)) {
                    sample_reports[sample_id] = []
                }
                sample_reports[sample_id].add(report)
            }
            // Emit per-sample tuples
            sample_reports.collect { sample_id, report_list ->
                def meta = meta_list.find { m -> m.id == sample_id }
                [meta, report_list]
            }
        }

    // Collect metawrap bins for batched GTDB-Tk
    ch_all_metawrap_bins = ch_metawrap.bins
        .collect()
        .map { items ->
            def meta_list = items.collect { item -> item[0] }
            def all_bins = []
            
            items.each { item ->
                def meta = item[0]
                def bins = item[1]
                all_bins.add(["${meta.id}", bins])
            }
            
            [meta_list, all_bins.collect { bin_tuple -> bin_tuple[1] }]
        }
    
    // Batched taxonomic classification with GTDB-TK
    ch_gtdb_tk_batch = GTDB_TK_BATCH(ch_all_metawrap_bins.combine(gtdbtk_db))
    ch_versions = ch_versions.mix(GTDB_TK_BATCH.out.versions.first())
    
    // Transform batched output back to per-sample format
    ch_gtdb_tk = ch_gtdb_tk_batch.report
        .flatMap { meta_list, reports ->
            // Group reports by sample
            def sample_reports = [:]
            reports.each { report ->
                def sample_id = report.name.split('_')[0]
                if (!sample_reports.containsKey(sample_id)) {
                    sample_reports[sample_id] = []
                }
                sample_reports[sample_id].add(report)
            }
            // Emit per-sample tuples
            sample_reports.collect { sample_id, report_list ->
                def meta = meta_list.find { m -> m.id == sample_id }
                [meta, report_list]
            }
        }

    //
    // Generate reports - mode-specific handling
    //
    if ( params.assembly_mode == "assembly" ) {
        // Per-sample mode: simple quality and taxonomy reports
        BIN_QUALITY_REPORT(ch_checkm_all_reports.map { _meta, reports -> reports }.collect())
        BIN_TAX_REPORT(ch_gtdb_tk.map { _meta, reports -> reports }.collect())
    } else if ( params.assembly_mode == "coassembly" ) {
        // Co-assembly mode: additional bin coverage analysis and summary report
        ch_bin_depth = BOWTIE2_SAMTOOLS_DEPTH(
            reads.combine(ch_metawrap.bins.map { _meta, bins -> bins })
        )
        ch_bin_cov = BEDTOOLS(ch_bin_depth)

        BIN_SUMMARY(
            ch_bin_cov.collect()
                .combine(ch_gtdb_tk.map { _meta, reports -> reports }.collect())
                .combine(ch_checkm_metawrap.map { _meta, reports -> reports }.collect())
        )
    }

    emit:
    refined_bins = ch_refined_bins    // channel: [ val(meta), path(bins) ] or path(bins)
    versions     = ch_versions        // channel: path(versions.yml)
}
