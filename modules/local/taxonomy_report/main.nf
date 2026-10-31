process TAXONOMY_REPORT {
    tag "$profiler-$db_name"
    label 'process_low'
    
    conda "conda-forge::python=3.13 conda-forge::pandas=3.1.0 conda-forge::numpy=2.5.4 conda-forge::matplotlib-base=3.11.2 conda-forge::seaborn-base=0.13.2 conda-forge::h5py=3.16.0 conda-forge::openpyxl=3.1.5 conda-forge::biopython=1.88"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:cd89ac4e30745df7' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:a1f6f67b03eef1d9' }"
    
    input:
    path(reports)
    // staged under a distinct name: the script's output is Reads_report.csv,
    // and writing through a same-named input symlink would corrupt the
    // upstream READS_REPORT work directory
    path(reads_report, stageAs: 'input_reads_report.csv')
    val(profiler)
    val(db_name)
    
    output:
    path("Reads_report.csv"), emit: report
    path("*.png"), emit: plots, optional: true
    path("versions.yml"), emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    taxonomy_report.py \\
        --profiler ${profiler} \\
        --reports ${reports} \\
        --reads-report ${reads_report} \\
        --db-name ${db_name} \\
        --output-dir . \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch Reads_report.csv
    touch versions.yml
    """
}
