process TAXONOMY_REPORT {
    tag "$profiler-$db_name"
    label 'process_low'
    
    conda "conda-forge::python=3.11 conda-forge::pandas=2.0 conda-forge::matplotlib=3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jupyter/scipy-notebook:python-3.11' :
        'jupyter/scipy-notebook:python-3.11' }"
    
    publishDir "${params.output}/reports/read_level/", mode: params.publish_dir_mode, pattern: 'Reads_report.csv'
    publishDir "${params.output}/reports/read_level/taxonomy", mode: params.publish_dir_mode, pattern: '*.png'
    
    input:
    path(reports)
    path(reads_report)
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
}
