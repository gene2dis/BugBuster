/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECKV_ENDTOEND Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Quality of geNomad's virus calls with CheckV end_to_end (roadmap Phase
    16, G2): completeness (from AAI to CheckV's genomes, else from its HMMs),
    host contamination, closed genomes, and the quality tier (Complete,
    High-quality, Medium-quality, Low-quality, Not-determined). The input is
    geNomad's <id>_filtered_contigs_virus.fna.gz: whole virus contigs and the
    provirus regions it found inside contigs (<contig>|provirus_<start>_<end>).
    CheckV only reads its database. A FASTA without sequences (geNomad found
    no virus) gives a header-only quality summary without running CheckV.

    Input:
        tuple val(meta), path(viruses) - geNomad's virus FASTA (gzipped)
        path(db)                       - CheckV database directory
                                         (genome_db/, hmm_db/)

    Output:
        results:         tuple val(meta), path("<id>_checkv") - quality_summary,
                         completeness, contamination, complete_genomes (and the
                         viruses.fna / proviruses.fna CheckV writes)
        quality_summary: tuple val(meta), path("<id>_checkv_quality_summary.tsv")
        versions:        path(versions.yml)
----------------------------------------------------------------------------------------
*/

process CHECKV_ENDTOEND {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.1.1--pyh106432d_1' :
        'quay.io/biocontainers/checkv:1.1.1--pyh106432d_1' }"

    input:
    tuple val(meta), path(viruses)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_checkv")                  , emit: results
    tuple val(meta), path("${meta.id}_checkv_quality_summary.tsv"), emit: quality_summary
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    """
    set -euo pipefail

    # count every header: 'grep -q' would stop at the first one and gzip die
    # of SIGPIPE, which pipefail reads as no sequences
    n_viruses=\$(gzip -cdf ${viruses} | grep -c '^>' || true)
    if [ "\$n_viruses" -gt 0 ]; then
        checkv end_to_end \\
            ${viruses} \\
            ${prefix}_checkv \\
            -d ${db} \\
            -t ${task.cpus} \\
            --remove_tmp \\
            ${args}
    else
        echo "No virus sequences in ${viruses}: header-only CheckV quality summary" >&2
        mkdir -p ${prefix}_checkv
        printf 'contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings\\n' \\
            > ${prefix}_checkv/quality_summary.tsv
    fi
    cp ${prefix}_checkv/quality_summary.tsv ${prefix}_checkv_quality_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv --version 2>&1 | sed 's/^CheckV, version //')
        diamond: \$(diamond --version 2>&1 | sed 's/^diamond version //')
        hmmer: \$(hmmsearch -h | sed -n 's/^# HMMER \\([0-9.]*\\).*/\\1/p')
        checkv_db: \$(sed -n 's/^## Ver[a-z]*on \\([0-9.]*\\).*/\\1/p' ${db}/README.txt | head -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p ${prefix}_checkv
    printf 'contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings\\n' \\
        > ${prefix}_checkv/quality_summary.tsv
    cp ${prefix}_checkv/quality_summary.tsv ${prefix}_checkv_quality_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: 1.1.1
        diamond: 2.1.8
        hmmer: 3.4
        checkv_db: 1.5
    END_VERSIONS
    """
}
