/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QUALITY CONTROL SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read quality filtering and host decontamination
----------------------------------------------------------------------------------------
*/

include { FASTP                      } from '../../modules/nf-core/fastp/main'
include { QFILTER                    } from '../../modules/local/qfilter/main'
include { COUNT_READS                } from '../../modules/local/count_reads/main'
include { BOWTIE2_DECONTAMINATE      } from '../../modules/local/bowtie2_decontaminate/main'
include { READS_REPORT               } from '../../modules/local/reads_report/main'

workflow QC {
    take:
    reads                   // channel: [ val(meta), [ reads ] ]
    decontamination_index   // channel: path(decontamination_index) - combined phiX + host index

    main:
    ch_versions = Channel.empty()
    ch_fastp_json = Channel.empty()
    
    if ( params.quality_control ) {
        //
        // Read quality filtering with nf-core Fastp
        // nf-core FASTP signature: tuple val(meta), path(reads), path(adapter_fasta)
        //                          val discard_trimmed_pass, val save_trimmed_fail, val save_merged
        //
        ch_fastp_input = reads.map { meta, reads_files -> 
            [ meta, reads_files, [] ]  // Empty adapter_fasta
        }
        
        FASTP(
            ch_fastp_input,
            false,  // discard_trimmed_pass
            false,  // save_trimmed_fail
            false   // save_merged
        )
        
        // Collect versions and FASTP json for MultiQC
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_fastp_json = FASTP.out.json
        
        //
        // Combine reads and json for QFILTER
        // QFILTER expects: tuple val(meta), path(reads), path(json)
        //
        ch_fastp_combined = FASTP.out.reads
            .join(FASTP.out.json)
            .map { meta, reads_files, json_file ->
                [ meta, reads_files, json_file ]
            }

        //
        // Extract and format QC reports
        //
        ch_fastp_reads_report = QFILTER(ch_fastp_combined)

        //
        // Filter samples by minimum read count
        //
        ch_fastp_reads_filtered = ch_fastp_reads_report.qfilter
            .map { meta, reads_out, read_count_file ->
                def after_reads = read_count_file.text.trim()
                if ( Integer.parseInt(after_reads) >= params.min_read_sample ) {
                    return [meta, reads_out]
                } else {
                    log.warn "Sample ${meta.id} has ${after_reads} reads after filtering, below threshold ${params.min_read_sample}. Skipping."
                    return null
                }
            }
            .filter { item -> item != null }

        //
        // Single-pass decontamination (removes phiX + host in one step)
        //
        ch_decontaminated = BOWTIE2_DECONTAMINATE(
            ch_fastp_reads_filtered.combine(decontamination_index),
            "contaminants"
        )

        //
        // Collect read reports
        //
        ch_reads_report = READS_REPORT(
            ch_decontaminated.report
                .concat(ch_fastp_reads_report.reads_report)
                .collect(),
            "contaminants"
        )

        ch_clean_reads = ch_decontaminated.reads
        ch_clean_reads_coassembly = ch_decontaminated.reads_coassembly
        ch_report = ch_reads_report

    } else {
        //
        // Skip QC - just count reads
        //
        ch_count = COUNT_READS(reads)
        
        ch_clean_reads = ch_count.reads
        ch_clean_reads_coassembly = ch_count.reads_coassembly
        ch_report = READS_REPORT(
            ch_count.reads_report.collect(),
            "none"
        )
    }

    emit:
    reads              = ch_clean_reads            // channel: [ val(meta), [ reads ] ]
    reads_coassembly   = ch_clean_reads_coassembly // channel: path(reads)
    report             = ch_report                 // channel: path(report)
    fastp_json         = ch_fastp_json             // channel: [ val(meta), path(json) ] for MultiQC
    versions           = ch_versions               // channel: path(versions.yml)
}
