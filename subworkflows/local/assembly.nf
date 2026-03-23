/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ASSEMBLY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Genome assembly using MEGAHIT with per-sample or co-assembly mode
    
    DESCRIPTION:
        This subworkflow handles metagenome assembly using MEGAHIT, followed by contig
        filtering with BBMap, and optional read alignment with Bowtie2/SAMtools for
        downstream binning and coverage analysis.
    
    INPUTS:
        reads: Channel of [meta, reads] tuples for per-sample assembly
        reads_coassembly: Channel of collected reads for co-assembly mode
    
    OUTPUTS:
        contigs: [meta, reads, contigs] - Filtered contigs with reads for binning
        contigs_meta: [meta, contigs] - Filtered contigs only for annotation
        bam: [meta, contigs, bam] - BAM files with contigs for depth analysis
        bam_meta: [meta, bam] - BAM files only for indexing
        versions: Tool version information
    
    MODES:
        assembly: Per-sample assembly (default)
        coassembly: Pool all reads and assemble together
        none: Skip assembly entirely
----------------------------------------------------------------------------------------
*/

include { MEGAHIT as NFCORE_MEGAHIT } from '../../modules/nf-core/megahit/main'
include { MEGAHIT              } from '../../modules/local/megahit/main'
include { BBMAP                } from '../../modules/local/bbmap/main'
include { BOWTIE2_SAMTOOLS     } from '../../modules/local/bowtie2_samtools/main'
include { CONTIG_FILTER_SUMMARY } from '../../modules/local/contig_filter_summary/main'

workflow ASSEMBLY {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    reads_coassembly   // channel: path(reads) - collected reads for coassembly

    main:
    ch_versions = Channel.empty()
    ch_contigs_with_reads = Channel.empty()
    ch_contigs_only = Channel.empty()
    ch_bam_with_contigs = Channel.empty()
    ch_bam_only = Channel.empty()
    ch_filter_reports = Channel.empty()
    ch_empty_reports = Channel.empty()

    //
    // Per-sample assembly mode
    //
    if ( params.assembly_mode == "assembly" ) {
        // Assemble reads with local MEGAHIT (unified process)
        MEGAHIT(reads)
        
        // Collect versions
        ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
        
        // Filter contigs by length with BBMap
        BBMAP(MEGAHIT.out.contigs_and_reads)
        
        // Collect BBMap versions
        ch_versions = ch_versions.mix(BBMAP.out.versions.first())
        
        // Collect filter reports for tracking
        ch_filter_reports = BBMAP.out.filter_report
        
        // Set output channels
        ch_contigs_with_reads = BBMAP.out.contigs_with_reads
        ch_contigs_only = BBMAP.out.contigs_only

        // Generate BAM files for binning/contig analysis if needed
        if ( params.include_binning || params.contig_tax_and_arg ) {
            BOWTIE2_SAMTOOLS(BBMAP.out.contigs_with_reads)
            
            ch_bam_with_contigs = BOWTIE2_SAMTOOLS.out.contigs_and_bam
            ch_bam_only = BOWTIE2_SAMTOOLS.out.bam_only
            ch_empty_reports = BOWTIE2_SAMTOOLS.out.empty_contig_report
            ch_versions = ch_versions.mix(BOWTIE2_SAMTOOLS.out.versions.first())
        }
    }

    //
    // Co-assembly mode
    //
    if ( params.assembly_mode == "coassembly" ) {
        // Prepare coassembly input: collect all reads from per-sample channel and pool them
        ch_coassembly_input = reads
            .map { _meta, reads_files -> reads_files }
            .collect()
            .map { all_reads -> [[id: "coassembly"], all_reads.flatten()] }
        
        // Co-assemble all reads with unified MEGAHIT process
        MEGAHIT(ch_coassembly_input)
        
        // Collect versions
        ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
        
        // Filter contigs by length with unified BBMAP process
        BBMAP(MEGAHIT.out.contigs_and_reads)
        
        // Collect BBMap versions
        ch_versions = ch_versions.mix(BBMAP.out.versions.first())
        
        // Collect filter reports for tracking
        ch_filter_reports = BBMAP.out.filter_report
        
        // Set output channels
        ch_contigs_only = BBMAP.out.contigs_only

        // Generate BAM files for binning/contig analysis if needed
        if ( params.include_binning || params.contig_tax_and_arg ) {
            // Prepare input: combine all original reads with filtered contigs
            // More efficient: collect reads first, then combine with contigs
            ch_alignment_input = reads
                .map { _meta, reads_files -> reads_files }
                .collect()
                .map { all_reads -> [[id: "coassembly"], all_reads.flatten()] }
                .combine(BBMAP.out.contigs_only.map { _meta, contigs -> contigs })
                .map { meta, reads_list, contigs -> [meta, reads_list, contigs] }
            
            BOWTIE2_SAMTOOLS(ch_alignment_input)
            
            ch_bam_only = BOWTIE2_SAMTOOLS.out.bam_only
            ch_bam_with_contigs = BOWTIE2_SAMTOOLS.out.contigs_and_bam
            ch_empty_reports = BOWTIE2_SAMTOOLS.out.empty_contig_report
            ch_versions = ch_versions.mix(BOWTIE2_SAMTOOLS.out.versions.first())
        }
    }

    // Generate contig filtering summary report
    if ( params.assembly_mode != "none" ) {
        CONTIG_FILTER_SUMMARY(
            ch_filter_reports.collect().ifEmpty([]),
            ch_empty_reports.collect().ifEmpty([])
        )
        ch_versions = ch_versions.mix(CONTIG_FILTER_SUMMARY.out.versions)
    }

    emit:
    contigs         = ch_contigs_with_reads  // channel: [ val(meta), path(reads), path(contigs) ] - for binning
    contigs_meta    = ch_contigs_only        // channel: [ val(meta), path(contigs) ] - for annotation
    bam             = ch_bam_with_contigs    // channel: [ val(meta), path(contigs), path(bam) ] - for depth analysis
    bam_meta        = ch_bam_only            // channel: [ val(meta), path(bam) ] - for indexing
    versions        = ch_versions            // channel: path(versions.yml)
}
