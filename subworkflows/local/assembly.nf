/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ASSEMBLY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Genome assembly using MEGAHIT with per-sample or co-assembly mode
----------------------------------------------------------------------------------------
*/

include { MEGAHIT as NFCORE_MEGAHIT } from '../../modules/nf-core/megahit/main'
include { MEGAHIT_COASSEMBLY   } from '../../modules/local/megahit/main'
include { BBMAP                } from '../../modules/local/bbmap/main'
include { BBMAP_COASSEMBLY     } from '../../modules/local/bbmap/main'
include { BOWTIE2_SAMTOOLS     } from '../../modules/local/bowtie2_samtools/main'
include { BOWTIE2_SAMTOOLS_COASSEMBLY } from '../../modules/local/bowtie2_samtools/main'

workflow ASSEMBLY {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    reads_coassembly   // channel: path(reads) - collected reads for coassembly

    main:
    ch_versions = Channel.empty()
    ch_contigs = Channel.empty()
    ch_contigs_meta = Channel.empty()
    ch_bam = Channel.empty()
    ch_bam_meta = Channel.empty()

    //
    // Per-sample assembly mode
    //
    if ( params.assembly_mode == "assembly" ) {
        //
        // Assemble reads with nf-core MEGAHIT
        // nf-core MEGAHIT signature:
        //   input:  tuple val(meta), path(reads1), path(reads2)
        //   output: tuple val(meta), path("*.contigs.fa.gz"), emit: contigs
        //
        ch_megahit_input = reads.map { meta, reads_files ->
            [ meta, reads_files[0], reads_files[1] ]  // Split reads into R1 and R2
        }
        
        NFCORE_MEGAHIT(ch_megahit_input)
        
        // Collect versions
        ch_versions = ch_versions.mix(NFCORE_MEGAHIT.out.versions.first())
        
        //
        // Re-join reads with contigs for BBMAP
        // Decompress contigs (nf-core outputs .fa.gz, BBMAP expects .fa)
        // Join original reads channel with contigs output
        //
        ch_assembly = reads
            .join(NFCORE_MEGAHIT.out.contigs)
            .map { meta, reads_files, contigs_gz ->
                [ meta, reads_files, contigs_gz ]
            }
        
        // Filter contigs by length with BBMap (handles gzipped input)
        ch_filtered = BBMAP(ch_assembly)
        
        ch_contigs = ch_filtered.bbmap
        ch_contigs_meta = ch_filtered.bbmap_contigs

        // Generate BAM files for binning/contig analysis
        if ( params.include_binning || params.contig_tax_and_arg ) {
            ch_bowtie2_samtools = BOWTIE2_SAMTOOLS(ch_filtered.bbmap)
            ch_bam = ch_bowtie2_samtools.contigs_and_bam
            ch_bam_meta = ch_bowtie2_samtools.only_bam
        }
    }

    //
    // Co-assembly mode
    //
    if ( params.assembly_mode == "coassembly" ) {
        // Co-assemble all reads with MEGAHIT
        ch_assembly_co = MEGAHIT_COASSEMBLY(reads_coassembly.collect())
        
        // Filter contigs by length
        ch_filtered_co = BBMAP_COASSEMBLY(ch_assembly_co)
        
        // Add meta information for coassembly
        ch_contigs_meta = ch_filtered_co.contigs.map { contigs ->
            def meta = [id: "coassembly"]
            return [meta, contigs]
        }

        // Generate BAM files for binning/contig analysis
        if ( params.include_binning || params.contig_tax_and_arg ) {
            ch_bowtie2_samtools_co = BOWTIE2_SAMTOOLS_COASSEMBLY(
                reads.combine(ch_filtered_co.contigs)
            )
            
            ch_bam_meta = ch_bowtie2_samtools_co.map { bam ->
                def meta = [id: "coassembly"]
                return [meta, bam]
            }
            
            ch_bam = ch_bowtie2_samtools_co
        }
    }

    emit:
    contigs         = ch_contigs          // channel: [ val(meta), path(reads), path(contigs) ]
    contigs_meta    = ch_contigs_meta     // channel: [ val(meta), path(contigs) ]
    bam             = ch_bam              // channel: path(bam) or [ val(meta), path(bam) ]
    bam_meta        = ch_bam_meta         // channel: [ val(meta), path(bam) ]
    versions        = ch_versions         // channel: path(versions.yml)
}
