/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BOWTIE2 Alignment Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read alignment and filtering using Bowtie2
----------------------------------------------------------------------------------------
*/

process BOWTIE2 {
    tag "${meta.id}_${db_alias}"
    container 'quay.io/biocontainers/bowtie2:2.5.3--py310ha0a81b8_0'

    label 'process_medium'

    publishDir "${params.output}/workflow/${meta.id}/clean_reads", 
        mode: params.publish_dir_mode,
        pattern: "*_clean_reads.fastq.gz"

    input:
    tuple val(meta), path(reads), path(index_db)
    val db_alias

    output:
    tuple val(meta), path("*map*fastq.gz"), emit: reads
    path("*_bowtie_report.tsv")           , emit: report
    path("*map*fastq.gz")                 , emit: reads_coassembly
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    def r1 = reads[0]
    def r2 = reads[1]
    def singleton = reads.size() > 2 ? reads[2] : null

    """
    # Get bowtie2 index basename
    basename=\$(ls ${index_db} | grep '\\.bt2' | sed -E 's/\\.rev\\.[0-9]\\.bt2//g' | sed -E 's/\\.[0-9]\\.bt2//g' | sort -u | head -1)

    # Align paired-end reads and extract unmapped
    bowtie2 \\
        ${args} \\
        -x ${index_db}/\${basename} \\
        -p ${task.cpus} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --un-conc-gz ${prefix}_${db_alias} > ${prefix}_bowtie_map_${db_alias}

    mv ${prefix}_${db_alias}.1 ${prefix}_R1_map_${db_alias}.fastq.gz
    mv ${prefix}_${db_alias}.2 ${prefix}_R2_map_${db_alias}.fastq.gz

    # Count reads
    R1_lines=\$(zcat ${prefix}_R1_map_${db_alias}.fastq.gz | wc -l)
    R1_count=\$(( R1_lines / 4 ))
    final_reads=\$(( R1_count * 2 ))

    # Process singleton reads if present
    ${singleton ? """
    bowtie2 \\
        ${args} \\
        -x ${index_db}/\${basename} \\
        -p ${task.cpus} \\
        -U ${singleton} \\
        --un-gz ${prefix}_singleton_unmapped > ${prefix}_bowtie_singleton_map_${db_alias}

    mv ${prefix}_singleton_unmapped ${prefix}_Singleton_map_${db_alias}.fastq.gz
    rm -f ${prefix}_bowtie_singleton_map_${db_alias}

    Singleton_lines=\$(zcat ${prefix}_Singleton_map_${db_alias}.fastq.gz | wc -l)
    Singleton_count=\$(( Singleton_lines / 4 ))
    echo -e \"Id\\tBowtie ${db_alias}\\tBowtie ${db_alias} singletons\" > ${prefix}_${db_alias}_bowtie_report.tsv
    echo -e \"${prefix}\\t\${final_reads}\\t\${Singleton_count}\" >> ${prefix}_${db_alias}_bowtie_report.tsv
    """ : """
    echo -e \"Id\\tBowtie ${db_alias}\" > ${prefix}_${db_alias}_bowtie_report.tsv
    echo -e \"${prefix}\\t\${final_reads}\" >> ${prefix}_${db_alias}_bowtie_report.tsv
    """}

    rm -f ${prefix}_bowtie_map_${db_alias}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch ${prefix}_R1_map_${db_alias}.fastq.gz
    touch ${prefix}_R2_map_${db_alias}.fastq.gz
    touch ${prefix}_${db_alias}_bowtie_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.3
    END_VERSIONS
    """
}
