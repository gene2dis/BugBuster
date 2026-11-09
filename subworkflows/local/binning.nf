/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BINNING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenomic binning using multiple tools and refinement with Binette
    
    DESCRIPTION:
        This subworkflow performs mode-agnostic metagenomic binning. It accepts
        contigs and BAM files (regardless of assembly or co-assembly origin) and
        runs a unified binning pipeline.
    
    INPUTS:
        contigs_and_bam: [meta, contigs, bam] - Contigs with aligned reads (BAM)
        checkm2_db: CheckM2 database for quality assessment
        gtdbtk_db: GTDB-Tk database for taxonomic classification
        sample_bams: [meta, bam] - One BAM per sample against its co-assembly
                (meta.group: the co-assembly's id; used only under co-assembly:
                MetaBAT2 depths and bin coverage)
        multisample_depth: [sample id, depth] - MetaBAT2's depth table with a
                column per sample of the group (MULTISAMPLE_COVERAGE, per-sample
                mode with --multisample_binning; empty otherwise)
    
    OUTPUTS:
        refined_bins: [meta, bins] - Binette-refined bins (the binner's bins with one binner)
        gtdbtk: [meta, reports] - The sample's GTDB-Tk bac120/ar53 summaries
        versions: Tool version information
    
    WORKFLOW:
        1. Calculate contig depth from BAM files (jgi; one column per sample
           under co-assembly): MetaBAT2's input and the MAG depth
        2. Run selected binning tools in parallel (from: MetaBAT2, SemiBin, COMEBin)
        3. If >=2 binners selected, refine bins with Binette
        4. Assess quality with CheckM2 (Binette's own run when refined)
        5. Classify taxonomy with GTDB-Tk
        6. MAG depth per sample (BIN_DEPTH, both modes) and reports
           (mode-specific: simple reports for assembly, bin summary for
           co-assembly)
----------------------------------------------------------------------------------------
*/

// Modules - unified for both per-sample and co-assembly modes
include { CALCULATE_DEPTH            } from '../../modules/local/calculate_depth/main'
include { METABAT2                   } from '../../modules/local/metabat2/main'
include { SEMIBIN                    } from '../../modules/local/semibin/main'
include { COMEBIN                    } from '../../modules/local/comebin/main'
include { BINETTE                    } from '../../modules/local/binette/main'
include { CHECKM2_BATCH              } from '../../modules/local/checkm2_batch/main'
include { GTDB_TK_BATCH              } from '../../modules/local/gtdb_tk_batch/main'
include { BIN_DEPTH                  } from '../../modules/local/bin_depth/main'
include { BIN_QUALITY_REPORT         } from '../../modules/local/bin_quality_report/main'
include { BIN_TAX_REPORT             } from '../../modules/local/bin_tax_report/main'
include { BIN_SUMMARY                } from '../../modules/local/bin_summary/main'
include { flagOn                     } from './utils_params'
include { multisampleSharedGroups    } from './utils_coassembly'

workflow BINNING {
    take:
    contigs_and_bam // channel: [ val(meta), path(contigs), path(bam) ]
    checkm2_db      // channel: path(checkm2_db)
    gtdbtk_db       // channel: path(gtdbtk_db)
    sample_bams     // channel: [ val(meta), path(bam) ] - per sample vs the co-assembly; used only under co-assembly
    multisample_depth // channel: [ val(sample id), path(depth) ] - --multisample_binning only

    main:
    ch_versions = channel.empty()

    // Parse selected binners
    def binners_list = params.binners instanceof List ? params.binners : params.binners.toString().tokenize(',').collect { b -> b.trim().toLowerCase() }
    def use_metabat2 = 'metabat2' in binners_list
    def use_semibin  = 'semibin' in binners_list
    def use_comebin  = 'comebin' in binners_list
    def num_binners  = binners_list.size()
    def use_refiner  = num_binners >= 2

    //
    // Unified binning workflow - mode-agnostic
    //
    
    // Contig depth from the BAM files with jgi_summarize_bam_contig_depths:
    // MetaBAT2's input and, with the refined bins, the MAG depth (BIN_DEPTH),
    // so it runs whatever the binners. Under co-assembly it gets one depth
    // column per sample, from the per-sample BAMs in sample order
    // (docs/benchmarks/binning_2026.md: MIMAG high-quality MetaBAT2 bins
    // 11 -> 23), and so does COMEBin (one coverage column per BAM:
    // docs/benchmarks/coassembly_comebin_bams_2026.md, roadmap Phase 15 B3);
    // SemiBin2 keeps the pooled BAM (its pre-trained models take one
    // coverage column)
    // Under co-assembly each co-assembly (one per --coassembly_group_by
    // group) takes its own samples' BAMs; it knows how many, so a group's
    // depths start once its BAMs are done
    ch_depth_input = params.assembly_mode == "coassembly"
        ? sample_bams
            .map { meta, bam -> [meta.group, meta.id, bam] }
            .combine(contigs_and_bam.map { meta, contigs, _bam -> [meta.id, meta, contigs] }, by: 0)
            .map { _label, id, bam, cmeta, contigs -> [groupKey([cmeta, contigs], cmeta.samples.size()), id, bam] }
            .groupTuple()
            .map { key, ids, bams ->
                def (cmeta, contigs) = key.getGroupTarget()
                [cmeta, contigs, [ids, bams].transpose().sort { a, b -> a[0] <=> b[0] }.collect { pair -> pair[1] }]
            }
        : contigs_and_bam
    CALCULATE_DEPTH(ch_depth_input)
    ch_versions = ch_versions.mix(CALCULATE_DEPTH.out.versions.first())
    // the samples each depth table has a column for: the co-assembly's, or the sample itself
    ch_depth_samples = CALCULATE_DEPTH.out.depth
        .map { meta, _contigs, depth -> [meta.id, depth, meta.coassembly ? meta.samples : [meta.id]] }

    // --multisample_binning (roadmap Phase 18b, B4): MetaBAT2 bins an
    // assembly of a group of two or more samples on the depth columns of every
    // sample of the group (MULTISAMPLE_COVERAGE); the MAG depth below keeps
    // the assembly's own table. Other assemblies bin as without the option
    ch_metabat_input = CALCULATE_DEPTH.out.depth
    if (flagOn(params.multisample_binning)) {
        def shared = multisampleSharedGroups()
        ch_metabat_input = CALCULATE_DEPTH.out.depth
            .filter { meta, _contigs, _depth -> !shared.containsKey(meta.id) }
            .mix(
                CALCULATE_DEPTH.out.depth
                    .filter { meta, _contigs, _depth -> shared.containsKey(meta.id) }
                    .map { meta, contigs, _depth -> [meta.id, meta, contigs] }
                    .join(multisample_depth)
                    .map { _id, meta, contigs, depth -> [meta, contigs, depth] }
            )
    }

    if (use_metabat2) {
        METABAT2(ch_metabat_input)
        ch_versions = ch_versions.mix(METABAT2.out.versions.first())
    }

    // Run selected binning tools
    if (use_semibin) {
        SEMIBIN(contigs_and_bam)
        ch_versions = ch_versions.mix(SEMIBIN.out.versions.first())
    }

    if (use_comebin) {
        COMEBIN(ch_depth_input)
        ch_versions = ch_versions.mix(COMEBIN.out.versions.first())
    }

    // Collect all selected binner outputs into a single channel
    ch_binner_bins = channel.empty()
    if (use_metabat2) ch_binner_bins = ch_binner_bins.mix(METABAT2.out.bins)
    if (use_semibin)  ch_binner_bins = ch_binner_bins.mix(SEMIBIN.out.bins)
    if (use_comebin)  ch_binner_bins = ch_binner_bins.mix(COMEBIN.out.bins)

    // Refine bins with Binette when >=2 binners are selected
    if (use_refiner) {
        // Group binner outputs by sample for Binette, with the sample's
        // contigs. groupTuple lists the bins in binner-completion order; sort
        // them by directory name so the task hash (and -resume) does not
        // depend on which binner finished first (design doc Q19)
        ch_bins_grouped = ch_binner_bins
            .map { meta, bins -> [meta.id, meta, bins] }
            .groupTuple(by: 0, size: num_binners)
            .map { id, metas, bins_list -> [id, metas[0], bins_list.sort(false) { b -> b.name }] }
            .join(contigs_and_bam.map { meta, contigs, _bam -> [meta.id, contigs] })
            .map { _id, meta, bins_list, contigs -> [meta, bins_list, contigs] }

        BINETTE(ch_bins_grouped.combine(checkm2_db))
        ch_refined_bins = BINETTE.out.bins
        ch_versions = ch_versions.mix(BINETTE.out.versions.first())

        // Binette scores every input and refined bin with CheckM2 and writes
        // the per-binner reports CHECKM2_BATCH would (performance review M7),
        // so the separate CheckM2 run is skipped
        ch_checkm_all_reports = BINETTE.out.quality_reports
        ch_checkm_refined = BINETTE.out.refined_report
    } else {
        // Single binner: use output directly, no refinement
        ch_refined_bins = ch_binner_bins

        // Collect all bins for batched CHECKM2 processing, sorted by sample and
        // bin directory so the batch input does not follow task completion order
        // (which changed the task hash and defeated -resume; design doc Q19)
        ch_all_sample_bins = ch_binner_bins
            .toSortedList { a, b -> (a[0].id <=> b[0].id) ?: (a[1].name <=> b[1].name) }
            .map { items ->
                def meta_list = items.collect { item -> item[0] }.unique { m -> m.id }
                def all_paths = items.collect { item -> item[1] }
                [meta_list, all_paths]
            }

        // Batched quality assessment with CheckM2
        ch_checkm_batch = CHECKM2_BATCH(ch_all_sample_bins.combine(checkm2_db))
        ch_versions = ch_versions.mix(CHECKM2_BATCH.out.versions.first())

        // Transform batched output back to per-sample format
        ch_checkm_all_reports = ch_checkm_batch.all_reports
            .flatMap { meta_list, reports ->
                // Group reports by sample
                def sample_reports = [:]
                def report_files = reports instanceof List ? reports : [reports]
                report_files.each { report ->
                    // Extract sample_id from filename: S10_Ago2021_metabat_quality_report.tsv -> S10_Ago2021
                    // Use findLastIndexOf to correctly handle sample names that contain binner keywords
                    def filename = report.name
                    def parts = filename.tokenize('_')
                    def binnerIdx = parts.findLastIndexOf { part -> part in ['metabat', 'semibin', 'comebin'] }
                    def sample_id = binnerIdx > 0 ? parts[0..<binnerIdx].join('_') : parts[0]

                    if (!sample_reports.containsKey(sample_id)) {
                        sample_reports[sample_id] = []
                    }
                    sample_reports[sample_id].add(report)
                }
                // Emit per-sample tuples
                sample_reports.collect { sample_id, report_list ->
                    def meta = meta_list.find { m ->
                        def mid = m instanceof Map ? (m.id ?: m['id']) : m.toString()
                        mid == sample_id
                    }
                    // Every report is named <sample>_<binner>_...: an unmatched
                    // one is a naming bug, and dropping it would silently skip
                    // the bin reports
                    if (meta == null) {
                        error "CheckM2: no sample matches the report(s) ${report_list*.name.join(', ')} (parsed sample id '${sample_id}')"
                    }
                    [meta, report_list]
                }
            }

        // The single binner's bins are the refined bins
        ch_checkm_refined = ch_checkm_all_reports
    }

    // Collect refined bins for batched GTDB-Tk (sorted, as for CHECKM2 above)
    ch_all_refined_bins = ch_refined_bins
        .filter { _meta, bins -> bins != null }
        .toSortedList { a, b -> (a[0].id <=> b[0].id) ?: (a[1].name <=> b[1].name) }
        .map { items ->
            def meta_list = items.collect { item -> item[0] }
            def all_paths = items.collect { item -> item[1] }
            [meta_list, all_paths]
        }
    
    // Batched taxonomic classification with GTDB-TK
    ch_gtdb_tk_batch = GTDB_TK_BATCH(ch_all_refined_bins.combine(gtdbtk_db))
    ch_versions = ch_versions.mix(GTDB_TK_BATCH.out.versions.first())
    
    // Transform batched output back to per-sample format
    ch_gtdb_tk = ch_gtdb_tk_batch.report
        .flatMap { meta_list, reports ->
            // Group reports by sample
            def sample_reports = [:]
            def report_files = reports instanceof List ? reports : [reports]
            report_files.each { report ->
                // Extract sample_id from filename: S10_Ago2021_gtdbtk_bac120.tsv -> S10_Ago2021
                // Use findLastIndexOf to correctly handle sample names that contain 'gtdbtk'
                def filename = report.name
                def parts = filename.tokenize('_')
                def gtdbIdx = parts.findLastIndexOf { part -> part == 'gtdbtk' }
                def sample_id = gtdbIdx > 0 ? parts[0..<gtdbIdx].join('_') : parts[0]
                
                if (!sample_reports.containsKey(sample_id)) {
                    sample_reports[sample_id] = []
                }
                sample_reports[sample_id].add(report)
            }
            // Emit per-sample tuples
            sample_reports.collect { sample_id, report_list ->
                def meta = meta_list.find { m -> 
                    def mid = m instanceof Map ? (m.id ?: m['id']) : m.toString()
                    mid == sample_id
                }
                if (meta == null) {
                    error "GTDB-Tk: no sample matches the report(s) ${report_list*.name.join(', ')} (parsed sample id '${sample_id}')"
                }
                [meta, report_list]
            }
        }

    //
    // MAG depth per sample, in both modes (roadmap Phase 15, B1): the
    // length-weighted mean of each refined bin's contig depths from the jgi
    // table MetaBAT2 bins on (alignments at >= 97 % identity). It needs at
    // least one refined bin: a bin set can be empty (Binette or SemiBin
    // write only SKIPPED.txt when nothing passes); BIN_SUMMARY then writes a
    // header-only table from the empty quality/taxonomy reports
    //
    def ch_bins_for_depth = ch_refined_bins
        .filter { meta, bins ->
            def has_bins = bins.listFiles().any { f -> f.name.endsWith('.fa') }
            if ( !has_bins ) {
                log.warn "No refined bins for ${meta.id}: skipping bin depth"
            }
            has_bins
        }
    BIN_DEPTH(
        ch_bins_for_depth
            .map { meta, bins -> [meta.id, meta, bins] }
            .join(ch_depth_samples)
            .map { _id, meta, bins, depth, samples -> [meta, depth, bins, samples] }
    )
    ch_versions = ch_versions.mix(BIN_DEPTH.out.versions.first())
    // <sample>_bin_depth.tsv, one per sample and assembly
    ch_bin_depth = BIN_DEPTH.out.bin_depth.map { _meta, tables -> tables }.flatten()

    //
    // Generate reports - mode-specific handling
    //
    if ( params.assembly_mode == "assembly" ) {
        // Per-sample mode: simple quality and taxonomy reports
        BIN_QUALITY_REPORT(ch_checkm_all_reports.map { _meta, reports -> reports }.collect(sort: true))
        BIN_TAX_REPORT(ch_gtdb_tk.map { _meta, reports -> reports }.collect(sort: true))
        ch_versions = ch_versions.mix(BIN_QUALITY_REPORT.out.versions, BIN_TAX_REPORT.out.versions)
    } else if ( params.assembly_mode == "coassembly" ) {
        // Co-assembly mode: the bin summary (quality, taxonomy and depth per
        // sample) of each co-assembly; without refined bins there are no depth
        // tables, and BIN_SUMMARY writes a header-only table from the empty
        // quality/taxonomy reports
        BIN_SUMMARY(
            ch_checkm_refined.map { meta, reports -> [meta.id, meta, reports] }
                .join(ch_gtdb_tk.map { meta, reports -> [meta.id, reports] })
                .map { id, meta, quality, taxonomy -> [id, [meta, [quality, taxonomy].flatten()]] }
                .join(BIN_DEPTH.out.bin_depth.map { meta, tables -> [meta.id, tables] }, remainder: true)
                .filter { _id, reports, _tables -> reports != null }
                .map { _id, reports, tables ->
                    def (meta, files) = reports
                    [meta, (files + [tables ?: []].flatten()).sort { f -> f.name }]
                }
        )
        ch_versions = ch_versions.mix(BIN_SUMMARY.out.versions)
    }

    emit:
    refined_bins = ch_refined_bins    // channel: [ val(meta), path(bins) ]
    gtdbtk       = ch_gtdb_tk         // channel: [ val(meta), [ GTDB-Tk bac120/ar53 summaries ] ]
    refined_quality = ch_checkm_refined  // channel: [ val(meta), CheckM2 report(s) of the refined bins ]
    bin_depth    = ch_bin_depth       // channel: path(<sample>_bin_depth.tsv), per sample and assembly
    versions     = ch_versions        // channel: path(versions.yml)
}
