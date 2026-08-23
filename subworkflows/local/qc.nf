/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QUALITY CONTROL SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read quality filtering and host decontamination
----------------------------------------------------------------------------------------
*/

include { FASTP                      } from '../../modules/nf-core/fastp/main'
include { FASTP as FASTP_SINGLETON   } from '../../modules/nf-core/fastp/main'
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
        // Paired reads and optional singletons are trimmed separately: the
        // nf-core FASTP paired branch only uses reads[0..1], so singletons get
        // their own single-end invocation and are rejoined afterwards.
        //
        ch_fastp_input = reads.map { meta, reads_files ->
            [ meta, reads_files[0..1], [] ]  // Empty adapter_fasta
        }

        ch_singleton_input = reads
            .filter { meta, reads_files -> reads_files.size() > 2 }
            .map { meta, reads_files ->
                [ meta + [single_end: true], [ reads_files[2] ], [] ]
            }

        FASTP(
            ch_fastp_input,
            false,  // discard_trimmed_pass
            false,  // save_trimmed_fail
            false   // save_merged
        )

        FASTP_SINGLETON(
            ch_singleton_input,
            false,  // discard_trimmed_pass
            false,  // save_trimmed_fail
            false   // save_merged
        )

        // Collect versions and FASTP json for MultiQC
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_versions = ch_versions.mix(FASTP_SINGLETON.out.versions.first())
        ch_fastp_json = FASTP.out.json.mix(FASTP_SINGLETON.out.json)

        //
        // Rejoin trimmed singletons with their paired reads, then combine with
        // the fastp jsons for QFILTER.
        // QFILTER expects: tuple val(meta), path(reads), path(json)
        //
        // The singleton side is packed into a single element so that
        // join(remainder: true) pads samples without singletons with exactly
        // one null.
        //
        ch_singleton_trimmed = FASTP_SINGLETON.out.reads
            .join(FASTP_SINGLETON.out.json)
            .map { meta, s_reads, s_json ->
                def s_read = s_reads instanceof List ? s_reads[0] : s_reads
                [ meta.id, [ s_read, s_json ] ]
            }

        ch_fastp_combined = FASTP.out.reads
            .join(FASTP.out.json)
            .map { meta, reads_files, json_file ->
                [ meta.id, meta, reads_files, json_file ]
            }
            .join(ch_singleton_trimmed, remainder: true)
            .map { _id, meta, reads_files, json_file, singleton ->
                singleton
                    ? [ meta, reads_files + [singleton[0]], [json_file, singleton[1]] ]
                    : [ meta, reads_files, json_file ]
            }

        //
        // Extract and format QC reports
        //
        ch_fastp_reads_report = QFILTER(ch_fastp_combined)

        //
        // Filter samples by minimum read count
        // (paired reads only: singletons do not count toward min_read_sample)
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
