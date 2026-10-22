/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BWAMEM3 DECONTAMINATE Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Single-pass read decontamination against a combined bwa-mem3 index
    (--read_aligner bwamem3). Removes reads mapping to any contaminant (phiX,
    host, etc.) in one step.

    bwa mem has no --un-conc/--un, so the clean reads are pulled from the
    alignment stream with samtools fastq. The pair policy is strict: a pair is
    clean only when BOTH mates are unmapped (-f 12); bowtie2 --un-conc kept
    any pair without a concordant alignment, so a pair with one host mate
    survived. -F 0x900 drops secondary/supplementary records, which would
    otherwise write the same read twice. -n keeps the input read names.
    -K fixes the batch size so the output does not depend on task.cpus.

    -P (no pairing) is required: in paired mode bwa mem applies the -T score
    floor only to single-end picks and emits a chosen pair without it. Random
    19-22 bp seed hits in a 3 Gb host then pass as "proper pairs" whenever the
    mates overlap (short inserts), which removed 17.5% of a bacterial mock
    community in the benchmark (0.1% with -P; simulated human reads 100%
    removed either way). Mate rescue stays on.
----------------------------------------------------------------------------------------
*/

process BWAMEM3_DECONTAMINATE {
    tag "${meta.id}"
    container 'community.wave.seqera.io/library/bwa-mem3_samtools:96d089551af9e9c3'

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
    set -euo pipefail

    # Index prefix: the *.bwt.2bit.64 file (bwa-mem3 and bwa-mem2 indexes)
    bwt=\$(find -L ${index_db} -maxdepth 1 -name '*.bwt.2bit.64' | sort | head -1)
    if [ -z "\${bwt}" ]; then
        echo "ERROR: no bwa-mem3/bwa-mem2 index (*.bwt.2bit.64) in ${index_db}." >&2
        if find -L ${index_db} -maxdepth 1 -name '*.bt2*' | grep -q .; then
            echo "It holds a Bowtie2 index: run with --read_aligner bowtie2 to keep using it, rebuild it with bwa-mem3, or pass the genome FASTA via --custom_host_fasta." >&2
        fi
        exit 1
    fi
    index=\${bwt%.bwt.2bit.64}

    # Align paired-end reads; keep pairs with both mates unmapped (clean).
    # -P: score each mate against -T on its own (see header)
    bwa-mem3 mem \\
        ${args} \\
        -P \\
        -t ${task.cpus} \\
        -K 100000000 \\
        --bam \\
        \${index} \\
        ${r1} \\
        ${r2} \\
        | samtools fastq \\
            -@ ${task.cpus} \\
            -n \\
            -f 12 \\
            -F 0x900 \\
            -1 ${prefix}_R1_clean.fastq.gz \\
            -2 ${prefix}_R2_clean.fastq.gz \\
            -0 /dev/null \\
            -s /dev/null \\
            - \\
            2> ${prefix}_fastq_paired.log \\
        || { rc=\$?; cat ${prefix}_fastq_paired.log >&2; exit \$rc; }
    cat ${prefix}_fastq_paired.log >&2

    # Clean reads = records that passed -f 12 -F 0x900 ("processed N reads")
    # minus those without their mate ("discarded N singletons", sent to -s),
    # i.e. the reads written to R1 + R2; no second pass over the output
    final_reads=\$(clean_read_count.sh ${prefix}_fastq_paired.log samtools)

    # Process singleton reads if present
    ${singleton ? """
    bwa-mem3 mem \\
        ${args} \\
        -t ${task.cpus} \\
        -K 100000000 \\
        --bam \\
        \${index} \\
        ${singleton} \\
        | samtools fastq \\
            -@ ${task.cpus} \\
            -n \\
            -f 4 \\
            -F 0x900 \\
            -0 ${prefix}_Singleton_clean.fastq.gz \\
            - \\
            2> ${prefix}_fastq_singleton.log \\
        || { rc=\$?; cat ${prefix}_fastq_singleton.log >&2; exit \$rc; }
    cat ${prefix}_fastq_singleton.log >&2

    Singleton_count=\$(clean_read_count.sh ${prefix}_fastq_singleton.log samtools)
    echo -e \"Id\\tClean reads (${db_alias} removed)\\tClean singletons (${db_alias} removed)\" > ${prefix}_${db_alias}_decontamination_report.tsv
    echo -e \"${prefix}\\t\${final_reads}\\t\${Singleton_count}\" >> ${prefix}_${db_alias}_decontamination_report.tsv
    """ : """
    echo -e \"Id\\tClean reads (${db_alias} removed)\" > ${prefix}_${db_alias}_decontamination_report.tsv
    echo -e \"${prefix}\\t\${final_reads}\" >> ${prefix}_${db_alias}_decontamination_report.tsv
    """}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem3: \$(bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
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
        bwamem3: 0.15.0
        samtools: 1.24
    END_VERSIONS
    """
}
