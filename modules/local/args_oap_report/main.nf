process ARGS_OAP_REPORT {
    tag "args_oap_report"
    label 'process_single'

    conda "conda-forge::python=3.13 conda-forge::pandas=3.1.0 conda-forge::numpy=2.5.4 conda-forge::matplotlib-base=3.11.2 conda-forge::seaborn-base=0.13.2 conda-forge::h5py=3.16.0 conda-forge::openpyxl=3.1.5 conda-forge::biopython=1.88 conda-forge::typst=0.15.1 conda-forge::jinja2=3.1.6 conda-forge::altair=6.3.0 conda-forge::vl-convert-python=1.9.0.post1 conda-forge::font-inter=4.1 conda-forge::pyyaml=6.0.3"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:2ce015cec686e756' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:6bff5ca407ba1454' }"

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
