process ARG_CONTIG_LEVEL_REPORT {
    tag "arg_contig"
    label 'process_single'
    
    conda "conda-forge::python=3.13 conda-forge::pandas=3.1.0 conda-forge::numpy=2.5.4 conda-forge::matplotlib-base=3.11.2 conda-forge::seaborn-base=0.13.2 conda-forge::h5py=3.16.0 conda-forge::openpyxl=3.1.5 conda-forge::biopython=1.88"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:cd89ac4e30745df7' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:a1f6f67b03eef1d9' }"

    input:
    path(arg_contig_data)

    output:
    path("Contig_tax_and_arg_prediction.tsv"), emit: arg_reports
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    arg_contig_level_report.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch Contig_tax_and_arg_prediction.tsv
    touch versions.yml
    """
}
