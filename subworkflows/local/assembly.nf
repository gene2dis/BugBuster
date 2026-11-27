/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ASSEMBLY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Genome assembly using MEGAHIT with per-sample or co-assembly mode
    
    DESCRIPTION:
        This subworkflow handles metagenome assembly using MEGAHIT, followed by contig
        filtering with BBMap, and optional read alignment with Bowtie2 or bwa-mem3
        (--read_aligner) plus SAMtools for downstream binning and coverage analysis.
        With --assembly_read_reduction normalize, BBNorm normalizes the reads
        MEGAHIT gets, and only those: BBMap, the alignments and everything
        downstream use the full reads.
    
    INPUTS:
        reads: Channel of [meta, reads] tuples for per-sample assembly
                (co-assembly mode pools this same channel internally: every
                sample into one co-assembly, or one co-assembly per
                --coassembly_group_by group)
    
    OUTPUTS (main.nf reads all of them):
        contigs_meta: [meta, contigs] - Filtered contigs (contig stages,
                multi-sample coverage)
        bam: [meta, contigs, bam] - Contigs with their BAM (BINNING's input)
        bam_meta: [meta, bam] - BAM files only (indexing, BlobTools)
        counting_bam: [meta, bam] - One BAM per sample, for gene quantification
                and co-assembly binning depths (assembly: the same BAM as
                bam_meta; coassembly: each sample against its co-assembly, whose
                samtools merge is bam/bam_meta - issue #28; meta.group names
                the co-assembly)
        stats: <id>_contig.stats - BBMap statistics of the MEGAHIT contigs
        filter_summary: contig_filtering_summary.txt
        empty_assemblies: one TSV line per assembly with no contig left
                after BBMAP's length filter; such an assembly is in none of
                the outputs above (roadmap Phase 21, RV10)
        versions: Tool version information

    CO-ASSEMBLY META:
        [id: 'coassembly' (or the group's label), coassembly: true,
         samples: [its sample ids, sorted]]
    
    MODES:
        assembly: Per-sample assembly (default)
        coassembly: Pool all reads and assemble together (per group with
                    --coassembly_group_by)
        none: Skip assembly entirely
----------------------------------------------------------------------------------------
*/

include { BBNORM               } from '../../modules/local/bbnorm/main'
include { MEGAHIT              } from '../../modules/local/megahit/main'
include { BBMAP                } from '../../modules/local/bbmap/main'
include { BOWTIE2_SAMTOOLS     } from '../../modules/local/bowtie2_samtools/main'
// --read_aligner bwamem3: the contig index is its own task, so the co-assembly
// index is built once and shared by the per-sample alignments
include { BWAMEM3_INDEX        } from '../../modules/local/bwamem3_samtools/main'
include { BWAMEM3_SAMTOOLS     } from '../../modules/local/bwamem3_samtools/main'
include { SAMTOOLS_MERGE as NFCORE_SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'
include { CONTIG_FILTER_SUMMARY } from '../../modules/local/contig_filter_summary/main'
include { flagOn; hasFastaRecord; bbmapMinLength } from './utils_params'
include { coassemblyGroups; coassemblyLabel } from './utils_coassembly'

workflow ASSEMBLY {
    take:
    reads              // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = channel.empty()
    ch_contigs_only = channel.empty()
    ch_bam_with_contigs = channel.empty()
    ch_bam_only = channel.empty()
    ch_sample_bams = channel.empty()
    ch_filter_reports = channel.empty()
    ch_contig_stats   = channel.empty()   // <id>_contig.stats (BBMap stats.sh)
    ch_filter_summary = channel.empty()   // contig_filtering_summary.txt
    ch_empty_reports = channel.empty()
    ch_empty_assemblies = channel.empty()  // metas of assemblies with no contig after BBMAP
    def reduce_reads = params.assembly_read_reduction == 'normalize'

    //
    // Per-sample assembly mode
    //
    if ( params.assembly_mode == "assembly" ) {
        // Assembler input: the clean reads, or their normalized copy
        ch_assembler_reads = reads
        if ( reduce_reads ) {
            BBNORM(reads)
            ch_assembler_reads = BBNORM.out.reads
            ch_versions = ch_versions.mix(BBNORM.out.versions.first())
        }

        // Assemble reads with local MEGAHIT (unified process)
        MEGAHIT(ch_assembler_reads)
        
        // Collect versions
        ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
        
        // Filter contigs by length with BBMap. BBMAP passes its reads on to
        // the alignment below, so with read reduction it gets the full reads
        // rather than MEGAHIT's input
        BBMAP(reduce_reads ? reads.join(MEGAHIT.out.contigs_only) : MEGAHIT.out.contigs_and_reads)
        
        // Collect BBMap versions
        ch_versions = ch_versions.mix(BBMAP.out.versions.first())
        
        // Collect filter reports for tracking
        ch_filter_reports = BBMAP.out.filter_report
        ch_contig_stats   = BBMAP.out.stats

        // An assembly with no contig left after the length filter goes no
        // further (roadmap Phase 21, RV10): see emptyAssemblies() below
        ch_kept = BBMAP.out.contigs_with_reads.branch { _meta, _reads, contigs ->
            kept: hasFastaRecord(contigs)
            empty: true
        }
        ch_empty_assemblies = ch_kept.empty.map { meta, _reads, _contigs -> meta }

        // Set output channels
        ch_contigs_with_reads = ch_kept.kept
        ch_contigs_only = ch_kept.kept.map { meta, _reads, contigs -> [meta, contigs] }

        // Generate BAM files for binning/contig analysis/gene counting if needed
        if ( flagOn(params.include_binning) || flagOn(params.contig_tax_and_arg) || flagOn(params.contig_level_functional) ) {
            if ( params.read_aligner == 'bwamem3' ) {
                BWAMEM3_INDEX(ch_contigs_only)
                ch_versions = ch_versions.mix(BWAMEM3_INDEX.out.versions.first())

                // [meta, reads, contigs] + [meta, index] -> [meta, reads, contigs, index]
                BWAMEM3_SAMTOOLS(ch_contigs_with_reads.join(BWAMEM3_INDEX.out.index))

                ch_bam_with_contigs = BWAMEM3_SAMTOOLS.out.contigs_and_bam
                ch_bam_only = BWAMEM3_SAMTOOLS.out.bam_only
                // Per-sample BAM is directly suitable for gene counting (Q1)
                ch_sample_bams = BWAMEM3_SAMTOOLS.out.bam_only
                ch_empty_reports = BWAMEM3_SAMTOOLS.out.empty_contig_report
                ch_versions = ch_versions.mix(BWAMEM3_SAMTOOLS.out.versions.first())
            } else {
                BOWTIE2_SAMTOOLS(ch_contigs_with_reads)

                ch_bam_with_contigs = BOWTIE2_SAMTOOLS.out.contigs_and_bam
                ch_bam_only = BOWTIE2_SAMTOOLS.out.bam_only
                // Per-sample BAM is directly suitable for gene counting (Q1)
                ch_sample_bams = BOWTIE2_SAMTOOLS.out.bam_only
                ch_empty_reports = BOWTIE2_SAMTOOLS.out.empty_contig_report
                ch_versions = ch_versions.mix(BOWTIE2_SAMTOOLS.out.versions.first())
            }
        }
    }

    //
    // Co-assembly mode
    //
    if ( params.assembly_mode == "coassembly" ) {
        // The co-assembly each sample goes into: its --coassembly_group_by
        // group's label, or 'coassembly' for every sample (roadmap Phase 15)
        def groups = coassemblyGroups()

        // Pool each co-assembly's reads, ordered by sample id (each sample
        // keeps its own R1, R2(, singleton) order, so MEGAHIT's R1/R2 lists
        // stay paired). groupTuple lists them in task completion order, which
        // changed the MEGAHIT task hash between identical runs and could re-run
        // the whole co-assembly on -resume (design doc Q19), hence the sort. It
        // has no group size: a sample QC drops (--min_read_sample) leaves its
        // group smaller, or empty, instead of waiting for it forever. Built once
        // for MEGAHIT (or BBNorm) and BBMAP.
        ch_coassembly_reads = reads
            .map { meta, files -> [coassemblyLabel(groups, meta.id), meta, files] }
            .groupTuple(by: 0)
            .map { label, metas, files_list ->
                def members = [metas, files_list].transpose().sort { a, b -> a[0].id <=> b[0].id }
                [
                    [id: label, coassembly: true, samples: members.collect { m -> m[0].id }],
                    members.collect { m -> m[1] }.flatten()
                ]
            }

        // Assembler input: the pooled reads, or their normalized copy (one
        // BBNorm call over every sample's reads)
        ch_assembler_reads = ch_coassembly_reads
        if ( reduce_reads ) {
            BBNORM(ch_coassembly_reads)
            ch_assembler_reads = BBNORM.out.reads
            ch_versions = ch_versions.mix(BBNORM.out.versions.first())
        }

        // Co-assemble all reads with unified MEGAHIT process
        MEGAHIT(ch_assembler_reads)
        
        // Collect versions
        ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
        
        // Filter contigs by length with unified BBMAP process (with the full
        // pooled reads under read reduction, as in per-sample mode)
        BBMAP(
            reduce_reads
                ? ch_coassembly_reads.join(MEGAHIT.out.contigs_only)
                : MEGAHIT.out.contigs_and_reads
        )
        
        // Collect BBMap versions
        ch_versions = ch_versions.mix(BBMAP.out.versions.first())
        
        // Collect filter reports for tracking
        ch_filter_reports = BBMAP.out.filter_report
        ch_contig_stats   = BBMAP.out.stats

        // A co-assembly with no contig left after the length filter goes no
        // further, and its samples get no alignment (roadmap Phase 21, RV10)
        ch_kept = BBMAP.out.contigs_only.branch { _meta, contigs ->
            kept: hasFastaRecord(contigs)
            empty: true
        }
        ch_empty_assemblies = ch_kept.empty.map { meta, _contigs -> meta }
        
        // Set output channels
        ch_contigs_only = ch_kept.kept

        // One alignment per sample against the co-assembly (issue #28): the
        // per-sample BAMs count genes (design doc Q1/T3), give MetaBAT2 and
        // COMEBin one coverage column per sample and bin depth per sample, and
        // their samtools merge is the pooled BAM SemiBin2 and BlobTools take. A merge
        // of the per-sample BAMs holds the same records as mapping every
        // sample in one task did. No read groups, as before: nothing
        // downstream reads them.
        if ( flagOn(params.include_binning) || flagOn(params.contig_tax_and_arg) || flagOn(params.contig_level_functional) ) {
            // [group label, sample meta + its group, reads, its co-assembly's contigs]
            ch_per_sample_input = reads
                .map { meta, files ->
                    def label = coassemblyLabel(groups, meta.id)
                    [label, meta + [group: label], files]
                }
                .combine(ch_contigs_only.map { meta, contigs -> [meta.id, contigs] }, by: 0)

            if ( params.read_aligner == 'bwamem3' ) {
                // One index per co-assembly, shared by its samples (bowtie2
                // builds it in each task)
                BWAMEM3_INDEX(ch_contigs_only)
                ch_versions = ch_versions.mix(BWAMEM3_INDEX.out.versions.first())

                BWAMEM3_SAMTOOLS(
                    ch_per_sample_input
                        .combine(BWAMEM3_INDEX.out.index.map { meta, index -> [meta.id, index] }, by: 0)
                        .map { _label, meta, files, contigs, index -> [meta, files, contigs, index] }
                )

                ch_sample_bams = BWAMEM3_SAMTOOLS.out.bam_only
                ch_empty_reports = BWAMEM3_SAMTOOLS.out.empty_contig_report
                ch_versions = ch_versions.mix(BWAMEM3_SAMTOOLS.out.versions.first())
            } else {
                BOWTIE2_SAMTOOLS(ch_per_sample_input.map { _label, meta, files, contigs -> [meta, files, contigs] })

                ch_sample_bams = BOWTIE2_SAMTOOLS.out.bam_only
                ch_empty_reports = BOWTIE2_SAMTOOLS.out.empty_contig_report
                ch_versions = ch_versions.mix(BOWTIE2_SAMTOOLS.out.versions.first())
            }
        }

        // Pooled BAM of each co-assembly for binning and contig taxonomy: its
        // samples' BAMs in sample order (a fixed task input, and the record
        // order of the old pooled task), named <id>_all_reads.bam (ext.prefix;
        // coassembly_all_reads.bam without groups). A co-assembly knows its
        // samples, so each group's merge starts once its own BAMs are done
        if ( flagOn(params.include_binning) || flagOn(params.contig_tax_and_arg) ) {
            NFCORE_SAMTOOLS_MERGE(
                ch_sample_bams
                    .map { meta, bam -> [meta.group, meta.id, bam] }
                    .combine(ch_contigs_only.map { meta, _contigs -> [meta.id, meta] }, by: 0)
                    .map { _label, id, bam, cmeta -> [groupKey(cmeta, cmeta.samples.size()), id, bam] }
                    .groupTuple()
                    .map { key, ids, bams ->
                        [key.getGroupTarget(), [ids, bams].transpose().sort { a, b -> a[0] <=> b[0] }.collect { pair -> pair[1] }, []]
                    },
                [[], [], [], []],
                ''
            )

            ch_bam_only = NFCORE_SAMTOOLS_MERGE.out.bam
            ch_bam_with_contigs = NFCORE_SAMTOOLS_MERGE.out.bam
                .join(ch_contigs_only)
                .map { meta, bam, contigs -> [meta, contigs, bam] }
        }
    }

    // Assemblies dropped for want of contigs (roadmap Phase 21, RV10; owner
    // decision 2026-11-21): like a sample below --min_read_sample, but only
    // for the contig stages. Its samples keep their QC, taxonomy and
    // read-level results; alignment, binning, contig taxonomy / ARGs and
    // contig functional annotation never see it (several of those tools stop
    // on an empty FASTA: MetaBAT2, BlobTools, DeepARG, eggNOG-mapper). One
    // line per assembly for pipeline_info/empty_assemblies.tsv: assembly,
    // its samples, contigs before the filter, the filter's minimum length
    ch_empty_rows = ch_empty_assemblies
        .map { meta -> [meta.id, meta] }
        .join(ch_filter_reports.map { report -> [report.name - '_filter_report.txt', report] })
        .map { id, meta, report ->
            def samples = meta.coassembly ? meta.samples : [meta.id]
            // both from BBMAP's own filter report: the minimum it applied may
            // come from ext.min_length rather than --bbmap_length
            def lines = report.readLines()
            def before = lines.find { line -> line.startsWith('contigs_before_filter\t') }?.tokenize('\t')?.getAt(1) ?: 'NA'
            def min_length = lines.find { line -> line.startsWith('min_length_threshold\t') }?.tokenize('\t')?.getAt(1) ?: bbmapMinLength()
            log.warn "Assembly ${id} (sample${samples.size() == 1 ? '' : 's'} ${samples.join(', ')}) has no contig of at least ${min_length} bp after the length filter (${before} contigs before it): it is dropped from alignment, binning, contig taxonomy and ARGs, and contig functional annotation; its samples keep their QC, taxonomy and read-level results (pipeline_info/empty_assemblies.tsv)"
            "${id}\t${samples.join(',')}\t${before}\t${min_length}\n"
        }

    // Generate contig filtering summary report
    if ( params.assembly_mode != "none" ) {
        CONTIG_FILTER_SUMMARY(
            ch_filter_reports.collect(sort: true).ifEmpty([]),
            ch_empty_reports.collect(sort: true).ifEmpty([])
        )
        ch_filter_summary = CONTIG_FILTER_SUMMARY.out.summary_report
        ch_versions = ch_versions.mix(CONTIG_FILTER_SUMMARY.out.versions)
    }

    emit:
    contigs_meta    = ch_contigs_only        // channel: [ val(meta), path(contigs) ] - for annotation
    bam             = ch_bam_with_contigs    // channel: [ val(meta), path(contigs), path(bam) ] - for depth analysis
    bam_meta        = ch_bam_only            // channel: [ val(meta), path(bam) ] - for indexing
    counting_bam    = ch_sample_bams         // channel: [ val(meta), path(bam) ] - per sample, against its assembly or the co-assembly
    stats           = ch_contig_stats        // channel: path(<id>_contig.stats) - BBMap statistics of the MEGAHIT contigs
    filter_summary  = ch_filter_summary      // channel: path(contig_filtering_summary.txt)
    empty_assemblies = ch_empty_rows         // channel: val("<assembly>\t<samples>\t<contigs before>\t<min length>\n"), assemblies dropped for want of contigs
    versions        = ch_versions            // channel: path(versions.yml)
}
