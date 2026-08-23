process TAXONOMY_REPORT {
    tag "$profiler-$db_name"
    label 'process_low'
    
    conda "conda-forge::python=3.11 conda-forge::pandas=2.0 conda-forge::numpy=1.24 conda-forge::matplotlib=3.7 conda-forge::seaborn=0.12 conda-forge::h5py=3.8 bioconda::biom-format=2.1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pruned:3b27e3935fded4e0' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pruned:195b3e3e5f741210' }"
    
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
