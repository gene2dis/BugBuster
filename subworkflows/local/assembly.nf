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
                (co-assembly mode pools this same channel internally)
    
    OUTPUTS:
        contigs: [meta, reads, contigs] - Filtered contigs with reads for binning
        contigs_meta: [meta, contigs] - Filtered contigs only for annotation
        bam: [meta, contigs, bam] - BAM files with contigs for depth analysis
        bam_meta: [meta, bam] - BAM files only for indexing
        counting_bam: [meta, bam] - One BAM per sample, for gene quantification
                and co-assembly binning depths (assembly: the same BAM as
                bam_meta; coassembly: each sample against the co-assembly, whose
                samtools merge is bam/bam_meta - issue #28)
        versions: Tool version information
    
    MODES:
        assembly: Per-sample assembly (default)
        coassembly: Pool all reads and assemble together
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
include { flagOn                } from './utils_params'

workflow ASSEMBLY {
    take:
    reads              // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = channel.empty()
    ch_contigs_with_reads = channel.empty()
    ch_contigs_only = channel.empty()
    ch_bam_with_contigs = channel.empty()
    ch_bam_only = channel.empty()
    ch_sample_bams = channel.empty()
    ch_filter_reports = channel.empty()
    ch_empty_reports = channel.empty()
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
        
        // Set output channels
        ch_contigs_with_reads = BBMAP.out.contigs_with_reads
        ch_contigs_only = BBMAP.out.contigs_only

        // Generate BAM files for binning/contig analysis/gene counting if needed
        if ( flagOn(params.include_binning) || flagOn(params.contig_tax_and_arg) || flagOn(params.contig_level_functional) ) {
            if ( params.read_aligner == 'bwamem3' ) {
                BWAMEM3_INDEX(BBMAP.out.contigs_only)
                ch_versions = ch_versions.mix(BWAMEM3_INDEX.out.versions.first())

                // [meta, reads, contigs] + [meta, index] -> [meta, reads, contigs, index]
                BWAMEM3_SAMTOOLS(BBMAP.out.contigs_with_reads.join(BWAMEM3_INDEX.out.index))

                ch_bam_with_contigs = BWAMEM3_SAMTOOLS.out.contigs_and_bam
                ch_bam_only = BWAMEM3_SAMTOOLS.out.bam_only
                // Per-sample BAM is directly suitable for gene counting (Q1)
                ch_sample_bams = BWAMEM3_SAMTOOLS.out.bam_only
                ch_empty_reports = BWAMEM3_SAMTOOLS.out.empty_contig_report
                ch_versions = ch_versions.mix(BWAMEM3_SAMTOOLS.out.versions.first())
            } else {
                BOWTIE2_SAMTOOLS(BBMAP.out.contigs_with_reads)

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
        // Pool all samples' reads for the co-assembly, ordered by sample id
        // (each sample keeps its own R1, R2(, singleton) order, so MEGAHIT's
        // R1/R2 lists stay paired). Plain collect() followed task completion
        // order, which changed the MEGAHIT task hash between identical runs
        // and could re-run the whole co-assembly on -resume (design doc Q19).
        // Built once for MEGAHIT (or BBNorm) and BBMAP.
        ch_pooled_reads = reads
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .map { items -> items.collect { item -> item[1] }.flatten() }

        ch_coassembly_reads = ch_pooled_reads.map { all_reads -> [[id: "coassembly"], all_reads] }

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
                ? ch_coassembly_reads.combine(MEGAHIT.out.contigs_only.map { _meta, contigs -> contigs })
                : MEGAHIT.out.contigs_and_reads
        )
        
        // Collect BBMap versions
        ch_versions = ch_versions.mix(BBMAP.out.versions.first())
        
        // Collect filter reports for tracking
        ch_filter_reports = BBMAP.out.filter_report
        
        // Set output channels
        ch_contigs_only = BBMAP.out.contigs_only

        // One alignment per sample against the co-assembly (issue #28): the
        // per-sample BAMs count genes (design doc Q1/T3), give MetaBAT2 one
        // depth column per sample and bin depth per sample, and their samtools
        // merge is the pooled BAM SemiBin2, COMEBin and BlobTools take. A merge
        // of the per-sample BAMs holds the same records as mapping every
        // sample in one task did. No read groups, as before: nothing
        // downstream reads them.
        if ( flagOn(params.include_binning) || flagOn(params.contig_tax_and_arg) || flagOn(params.contig_level_functional) ) {
            ch_per_sample_input = reads
                .combine(BBMAP.out.contigs_only.map { _meta, contigs -> contigs })
                .map { meta, reads_files, contigs -> [meta, reads_files, contigs] }

            if ( params.read_aligner == 'bwamem3' ) {
                // One index of the co-assembly, shared by every sample (bowtie2
                // builds it in each task)
                BWAMEM3_INDEX(BBMAP.out.contigs_only)
                ch_versions = ch_versions.mix(BWAMEM3_INDEX.out.versions.first())

                BWAMEM3_SAMTOOLS(ch_per_sample_input.combine(BWAMEM3_INDEX.out.index.map { _meta, index -> index }))

                ch_sample_bams = BWAMEM3_SAMTOOLS.out.bam_only
                ch_empty_reports = BWAMEM3_SAMTOOLS.out.empty_contig_report
                ch_versions = ch_versions.mix(BWAMEM3_SAMTOOLS.out.versions.first())
            } else {
                BOWTIE2_SAMTOOLS(ch_per_sample_input)

                ch_sample_bams = BOWTIE2_SAMTOOLS.out.bam_only
                ch_empty_reports = BOWTIE2_SAMTOOLS.out.empty_contig_report
                ch_versions = ch_versions.mix(BOWTIE2_SAMTOOLS.out.versions.first())
            }
        }

        // Pooled BAM for binning and contig taxonomy: the per-sample BAMs in
        // sample order (a fixed task input, and the record order of the old
        // pooled task), named coassembly_all_reads.bam (ext.prefix)
        if ( flagOn(params.include_binning) || flagOn(params.contig_tax_and_arg) ) {
            NFCORE_SAMTOOLS_MERGE(
                ch_sample_bams
                    .toSortedList { a, b -> a[0].id <=> b[0].id }
                    .map { items -> [[id: "coassembly"], items.collect { item -> item[1] }, []] },
                [[], [], [], []],
                ''
            )

            ch_bam_only = NFCORE_SAMTOOLS_MERGE.out.bam
            ch_bam_with_contigs = NFCORE_SAMTOOLS_MERGE.out.bam
                .combine(BBMAP.out.contigs_only.map { _meta, contigs -> contigs })
                .map { meta, bam, contigs -> [meta, contigs, bam] }
        }
    }

    // Generate contig filtering summary report
    if ( params.assembly_mode != "none" ) {
        CONTIG_FILTER_SUMMARY(
            ch_filter_reports.collect(sort: true).ifEmpty([]),
            ch_empty_reports.collect(sort: true).ifEmpty([])
        )
        ch_versions = ch_versions.mix(CONTIG_FILTER_SUMMARY.out.versions)
    }

    emit:
    contigs         = ch_contigs_with_reads  // channel: [ val(meta), path(reads), path(contigs) ] - for binning
    contigs_meta    = ch_contigs_only        // channel: [ val(meta), path(contigs) ] - for annotation
    bam             = ch_bam_with_contigs    // channel: [ val(meta), path(contigs), path(bam) ] - for depth analysis
    bam_meta        = ch_bam_only            // channel: [ val(meta), path(bam) ] - for indexing
    counting_bam    = ch_sample_bams         // channel: [ val(meta), path(bam) ] - per sample, against its assembly or the co-assembly
    versions        = ch_versions            // channel: path(versions.yml)
}
