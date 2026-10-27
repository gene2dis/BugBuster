process ARGS_OAP_REPORT {
    tag "args_oap_report"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    path(args_oap_dirs)

    output:
    path("args_oap_*.tsv"), emit: reports
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    args_oap_report.py ${args_oap_dirs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    printf 'sample\\tlevel\\tsarg_type\\tsarg_subtype\\tsarg_gene\\tread_count\\tcopies\\tcopies_per_cell\\tcopies_per_16S\\tppm\\n' > args_oap_long.tsv
    printf 'sample\\tnRead\\tn16S\\tnCell\\tsarg_read_count\\tstatus\\n' > args_oap_cells.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 1.5.3
    END_VERSIONS
    """
}
