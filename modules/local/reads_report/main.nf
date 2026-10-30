process READS_REPORT {

    conda "conda-forge::python=3.13 conda-forge::pandas=3.1.0 conda-forge::numpy=2.5.4 conda-forge::matplotlib-base=3.11.2 conda-forge::seaborn-base=0.13.2 conda-forge::h5py=3.16.0 conda-forge::openpyxl=3.1.5 conda-forge::biopython=1.88"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:cd89ac4e30745df7' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:a1f6f67b03eef1d9' }"

    tag "reads_report"
    label 'process_single'

    input:
        path(reports)
        val(args)

    output:
        path("Reads_report.csv"), emit: report
        path("Box_plot_reads.png"), emit: plot, optional: true
        path "versions.yml", emit: versions

    script:
        """
        report_unify.py ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed 's/Python //')
            pandas: \$(python3 -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo unknown)
        END_VERSIONS
	"""

    stub:
        """
        touch Reads_report.csv
        touch Box_plot_reads.png

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.13.18
            pandas: 3.1.0
        END_VERSIONS
        """
}
