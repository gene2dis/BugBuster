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

    input:
    tuple val(meta), path(reads), path(index_db)
    val db_alias

    output:
    // Multi-file glob outputs are sorted lexicographically: [R1, R2, Singleton]
    tuple val(meta), path("${meta.id}_{R1,R2,Singleton}_clean.fastq.gz"), emit: reads
    path("*_decontamination_report.tsv")                                , emit: report
    path("${meta.id}_*_clean.fastq.gz")                                 , emit: reads_coassembly
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    def r1 = reads[0]
    def r2 = reads[1]
    def singleton = reads.size() > 2 ? reads[2] : null

    """
    # Get bowtie2 index basename (.bt2, or .bt2l for a large index: bowtie2-build
    # writes one for references over ~4 Gbp, e.g. several host genomes)
    basename=\$(ls ${index_db} | grep -E '\\.[0-9]\\.bt2l?\$' | sed -E 's/(\\.rev)?\\.[0-9]\\.bt2l?\$//' | sort -u | head -1 || true)
    if [ -z "\${basename}" ]; then
        echo "ERROR: no Bowtie2 index (*.bt2 / *.bt2l) in ${index_db}." >&2
        if find -L ${index_db} -maxdepth 1 -name '*.bwt.2bit.64' | grep -q .; then
            echo "It holds a bwa-mem3/bwa-mem2 index: run with the default --read_aligner bwamem3 to use it." >&2
        fi
        exit 1
    fi

    # Align paired-end reads and extract unmapped (clean) reads
    bowtie2 \\
        ${args} \\
        -x ${index_db}/\${basename} \\
        -p ${task.cpus} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --un-conc-gz ${prefix}_clean_tmp.%.fastq.gz > ${prefix}_bowtie_mapped_${db_alias}

    # Rename output files. bowtie2 puts the mate number at the '%'; without
    # one it inserts it before the last extension, which for a dotted sample
    # id (S1.2) lands inside the id (S1.1.2_clean_tmp)
    mv ${prefix}_clean_tmp.1.fastq.gz ${prefix}_R1_clean.fastq.gz
    mv ${prefix}_clean_tmp.2.fastq.gz ${prefix}_R2_clean.fastq.gz

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
    def has_singleton = reads.size() > 2
    """
    echo -n '' | gzip > ${prefix}_R1_clean.fastq.gz
    echo -n '' | gzip > ${prefix}_R2_clean.fastq.gz
    ${has_singleton ? "echo -n '' | gzip > ${prefix}_Singleton_clean.fastq.gz" : ""}

    ${has_singleton ? """
    echo -e \"Id\\tClean reads (${db_alias} removed)\\tClean singletons (${db_alias} removed)\" > ${prefix}_${db_alias}_decontamination_report.tsv
    echo -e \"${prefix}\\t1000\\t100\" >> ${prefix}_${db_alias}_decontamination_report.tsv
    """ : """
    echo -e \"Id\\tClean reads (${db_alias} removed)\" > ${prefix}_${db_alias}_decontamination_report.tsv
    echo -e \"${prefix}\\t1000\" >> ${prefix}_${db_alias}_decontamination_report.tsv
    """}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.3
    END_VERSIONS
    """
}
