/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QUALITY CONTROL SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read quality filtering and host decontamination
----------------------------------------------------------------------------------------
*/

include { FASTP                      } from '../../modules/fastp/main'
include { QFILTER                    } from '../../modules/qfilter/main'
include { COUNT_READS                } from '../../modules/count_reads/main'
include { BOWTIE2 as BOWTIE2_PHIX    } from '../../modules/bowtie2/main'
include { BOWTIE2 as BOWTIE2_HOST    } from '../../modules/bowtie2/main'
include { READS_REPORT               } from '../../modules/reads_report/main'

workflow QC {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    phix_index      // channel: path(phix_index)
    host_index      // channel: path(host_index)

    main:
    ch_versions = Channel.empty()
    
    if ( params.quality_control ) {
        //
        // Read quality filtering with Fastp
        //
        ch_fastp_reads = FASTP(reads)

        //
        // Extract and format QC reports
        //
        ch_fastp_reads_report = QFILTER(ch_fastp_reads)

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
            .filter { it != null }

        //
        // Remove PhiX contamination
        //
        ch_phix_clean = BOWTIE2_PHIX(
            ch_fastp_reads_filtered.combine(phix_index),
            "phiX"
        )

        //
        // Remove host contamination
        //
        ch_host_clean = BOWTIE2_HOST(
            ch_phix_clean.reads.combine(host_index),
            "host"
        )

        //
        // Collect read reports
        //
        ch_reads_report = READS_REPORT(
            ch_phix_clean.report
                .concat(ch_host_clean.report)
                .concat(ch_fastp_reads_report.reads_report)
                .collect(),
            "phiX host"
        )

        ch_clean_reads = ch_host_clean.reads
        ch_clean_reads_coassembly = ch_host_clean.reads_coassembly
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
    versions           = ch_versions               // channel: path(versions.yml)
}
