process KMARVEL_REPORT {
    tag "kmarvel_report"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:cd89ac4e30745df7' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:a1f6f67b03eef1d9' }"

    input:
    path(kmarvel_tsvs)

    output:
    path("KMARVEL_*.tsv"), emit: reports
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    kmarvel_report.py ${kmarvel_tsvs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch KMARVEL_summary_long.tsv KMARVEL_supporting_reads_matrix.tsv KMARVEL_gene_family_matrix.tsv \\
        KMARVEL_drug_class_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 1.5.3
    END_VERSIONS
    """
}
