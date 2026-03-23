/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BOWTIE2 DECONTAMINATE Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Single-pass read decontamination using combined Bowtie2 index
    Removes reads mapping to any contaminant (phiX, host, etc.) in one step
----------------------------------------------------------------------------------------
*/

process BOWTIE2_DECONTAMINATE {
    tag "${meta.id}"
    container 'quay.io/biocontainers/bowtie2:2.5.3--py310ha0a81b8_0'

    label 'process_high'

    storeDir params.store_clean_reads ? "${params.output}/clean_reads/${meta.id}" : null

    input:
    tuple val(meta), path(reads), path(index_db)
    val db_alias

    output:
    tuple val(meta), path("${meta.id}_R*_clean.fastq.gz"), emit: reads
    path("*_decontamination_report.tsv")                  , emit: report
    path("${meta.id}_*_clean.fastq.gz")                   , emit: reads_coassembly
    path "versions.yml"                                   , emit: versions

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

    # Align paired-end reads and extract unmapped (clean) reads
    bowtie2 \\
        ${args} \\
        -x ${index_db}/\${basename} \\
        -p ${task.cpus} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --un-conc-gz ${prefix}_clean_tmp > ${prefix}_bowtie_mapped_${db_alias}

    # Rename output files
    mv ${prefix}_clean_tmp.1 ${prefix}_R1_clean.fastq.gz
    mv ${prefix}_clean_tmp.2 ${prefix}_R2_clean.fastq.gz

    # Count clean reads
    R1_lines=\$(zcat ${prefix}_R1_clean.fastq.gz | wc -l)
    R1_count=\$(( R1_lines / 4 ))
    final_reads=\$(( R1_count * 2 ))

    # Process singleton reads if present
    ${singleton ? """
    bowtie2 \\
        ${args} \\
        -x ${index_db}/\${basename} \\
        -p ${task.cpus} \\
        -U ${singleton} \\
        --un-gz ${prefix}_singleton_clean > ${prefix}_bowtie_singleton_mapped_${db_alias}

    mv ${prefix}_singleton_clean ${prefix}_Singleton_clean.fastq.gz
    rm -f ${prefix}_bowtie_singleton_mapped_${db_alias}

    Singleton_lines=\$(zcat ${prefix}_Singleton_clean.fastq.gz | wc -l)
    Singleton_count=\$(( Singleton_lines / 4 ))
    echo -e \"Id\\tClean reads (${db_alias} removed)\\tClean singletons (${db_alias} removed)\" > ${prefix}_${db_alias}_decontamination_report.tsv
    echo -e \"${prefix}\\t\${final_reads}\\t\${Singleton_count}\" >> ${prefix}_${db_alias}_decontamination_report.tsv
    """ : """
    echo -e \"Id\\tClean reads (${db_alias} removed)\" > ${prefix}_${db_alias}_decontamination_report.tsv
    echo -e \"${prefix}\\t\${final_reads}\" >> ${prefix}_${db_alias}_decontamination_report.tsv
    """}

    # Clean up alignment output
    rm -f ${prefix}_bowtie_mapped_${db_alias}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch ${prefix}_R1_clean.fastq.gz
    touch ${prefix}_R2_clean.fastq.gz
    touch ${prefix}_${db_alias}_decontamination_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.3
    END_VERSIONS
    """
}
