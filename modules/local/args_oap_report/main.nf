process ARGS_OAP_REPORT {
    tag "args_oap_report"
    label 'process_single'

    conda "conda-forge::python=3.13 conda-forge::pandas=3.1.0 conda-forge::numpy=2.5.4 conda-forge::matplotlib-base=3.11.2 conda-forge::seaborn-base=0.13.2 conda-forge::h5py=3.16.0 conda-forge::openpyxl=3.1.5 conda-forge::biopython=1.88"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:cd89ac4e30745df7' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:a1f6f67b03eef1d9' }"

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
        python: 3.13.18
        pandas: 3.1.0
    END_VERSIONS
    """
}
