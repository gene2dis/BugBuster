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
        counting_bam: [meta, bam] - Per-sample BAMs for gene quantification
                (assembly: the same per-sample BAM as bam_meta; coassembly:
                dedicated per-sample alignments against the co-assembly, since
                the pooled BAM has no read groups - design doc Q1/T3)
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
// Aliased second instance: per-sample alignments against the co-assembly for
// gene quantification (design doc Q1/T3). Per-sample meta.id selects the
// module's per-sample script branch; the coassembly index is rebuilt in each
// task (accepted trade-off, revisit at T9 if it dominates runtime).
include { BOWTIE2_SAMTOOLS as BOWTIE2_SAMTOOLS_PER_SAMPLE } from '../../modules/local/bowtie2_samtools/main'
// --read_aligner bwamem3: the contig index is its own task, so the co-assembly
// index is built once and shared by the pooled and per-sample alignments
include { BWAMEM3_INDEX        } from '../../modules/local/bwamem3_samtools/main'
include { BWAMEM3_SAMTOOLS     } from '../../modules/local/bwamem3_samtools/main'
include { BWAMEM3_SAMTOOLS as BWAMEM3_SAMTOOLS_PER_SAMPLE } from '../../modules/local/bwamem3_samtools/main'
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
    ch_counting_bam = channel.empty()
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
                ch_counting_bam = BWAMEM3_SAMTOOLS.out.bam_only
                ch_empty_reports = BWAMEM3_SAMTOOLS.out.empty_contig_report
                ch_versions = ch_versions.mix(BWAMEM3_SAMTOOLS.out.versions.first())
            } else {
                BOWTIE2_SAMTOOLS(BBMAP.out.contigs_with_reads)

                ch_bam_with_contigs = BOWTIE2_SAMTOOLS.out.contigs_and_bam
                ch_bam_only = BOWTIE2_SAMTOOLS.out.bam_only
                // Per-sample BAM is directly suitable for gene counting (Q1)
                ch_counting_bam = BOWTIE2_SAMTOOLS.out.bam_only
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
        // Built once and reused for the co-assembly alignment below.
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

        // One bwa-mem3 index of the co-assembly, shared by the pooled and the
        // per-sample alignments below (bowtie2 rebuilds it in each task)
        ch_coassembly_index = channel.empty()
        if ( params.read_aligner == 'bwamem3' && ( flagOn(params.include_binning) || flagOn(params.contig_tax_and_arg) || flagOn(params.contig_level_functional) ) ) {
            BWAMEM3_INDEX(BBMAP.out.contigs_only)
            ch_coassembly_index = BWAMEM3_INDEX.out.index.map { _meta, index -> index }
            ch_versions = ch_versions.mix(BWAMEM3_INDEX.out.versions.first())
        }

        // Generate BAM files for binning/contig analysis if needed
        if ( flagOn(params.include_binning) || flagOn(params.contig_tax_and_arg) ) {
            // Prepare input: combine all original reads with filtered contigs
            // More efficient: collect reads first, then combine with contigs
            ch_alignment_input = ch_coassembly_reads
                .combine(BBMAP.out.contigs_only.map { _meta, contigs -> contigs })
                .map { meta, reads_list, contigs -> [meta, reads_list, contigs] }

            if ( params.read_aligner == 'bwamem3' ) {
                BWAMEM3_SAMTOOLS(ch_alignment_input.combine(ch_coassembly_index))

                ch_bam_only = BWAMEM3_SAMTOOLS.out.bam_only
                ch_bam_with_contigs = BWAMEM3_SAMTOOLS.out.contigs_and_bam
                ch_empty_reports = BWAMEM3_SAMTOOLS.out.empty_contig_report
                ch_versions = ch_versions.mix(BWAMEM3_SAMTOOLS.out.versions.first())
            } else {
                BOWTIE2_SAMTOOLS(ch_alignment_input)

                ch_bam_only = BOWTIE2_SAMTOOLS.out.bam_only
                ch_bam_with_contigs = BOWTIE2_SAMTOOLS.out.contigs_and_bam
                ch_empty_reports = BOWTIE2_SAMTOOLS.out.empty_contig_report
                ch_versions = ch_versions.mix(BOWTIE2_SAMTOOLS.out.versions.first())
            }
        }

        // Per-sample alignments against the co-assembly for gene counting
        // (design doc Q1/T3): the pooled BAM above has no read groups, so it
        // cannot yield per-sample counts. Per-sample meta.id (!= 'coassembly')
        // selects the module's per-sample script branch.
        if ( flagOn(params.contig_level_functional) ) {
            ch_per_sample_input = reads
                .combine(BBMAP.out.contigs_only.map { _meta, contigs -> contigs })
                .map { meta, reads_files, contigs -> [meta, reads_files, contigs] }

            if ( params.read_aligner == 'bwamem3' ) {
                BWAMEM3_SAMTOOLS_PER_SAMPLE(ch_per_sample_input.combine(ch_coassembly_index))

                ch_counting_bam = BWAMEM3_SAMTOOLS_PER_SAMPLE.out.bam_only
                ch_empty_reports = ch_empty_reports.mix(BWAMEM3_SAMTOOLS_PER_SAMPLE.out.empty_contig_report)
                ch_versions = ch_versions.mix(BWAMEM3_SAMTOOLS_PER_SAMPLE.out.versions.first())
            } else {
                BOWTIE2_SAMTOOLS_PER_SAMPLE(ch_per_sample_input)

                ch_counting_bam = BOWTIE2_SAMTOOLS_PER_SAMPLE.out.bam_only
                ch_empty_reports = ch_empty_reports.mix(BOWTIE2_SAMTOOLS_PER_SAMPLE.out.empty_contig_report)
                ch_versions = ch_versions.mix(BOWTIE2_SAMTOOLS_PER_SAMPLE.out.versions.first())
            }
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
    counting_bam    = ch_counting_bam        // channel: [ val(meta), path(bam) ] - per-sample, for gene quantification
    versions        = ch_versions            // channel: path(versions.yml)
}
