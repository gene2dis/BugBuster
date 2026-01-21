/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTP Quality Filtering Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read quality filtering and adapter trimming
----------------------------------------------------------------------------------------
*/

process FASTP {
    tag "${meta.id}"
    container 'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0'

    label 'process_low'

    publishDir "${params.output}/workflow/${meta.id}/fastp", mode: params.publish_dir_mode, pattern: "*_report.json"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastp.fastq.gz"), path("*_report.json"), emit: fastq
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    def r1 = reads[0]
    def r2 = reads[1]
    def singleton = reads.size() > 2 ? reads[2] : null

    """
    # Calculate optimal thread usage
    threads=\$(( ${task.cpus} > 1 ? ${task.cpus} / 2 : 1 ))

    # Process paired-end reads
    fastp \\
        --in1 ${r1} \\
        --in2 ${r2} \\
        --out1 ${prefix}_R1_fastp.fastq.gz \\
        --out2 ${prefix}_R2_fastp.fastq.gz \\
        -j ${prefix}_report.json \\
        ${args} \\
        --thread \$threads

    # Process singleton reads if present
    ${singleton ? """
    fastp \\
        -i ${singleton} \\
        -o ${prefix}_Singleton_fastp.fastq.gz \\
        -j ${prefix}_Singleton_report.json \\
        ${args} \\
        --thread \$threads
    """ : ""}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch ${prefix}_R1_fastp.fastq.gz
    touch ${prefix}_R2_fastp.fastq.gz
    touch ${prefix}_report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: 0.23.2
    END_VERSIONS
    """
}
