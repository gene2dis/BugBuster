process READS_REPORT {

    container 'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0'

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
            python: 3.9.12
            pandas: 1.4.2
        END_VERSIONS
        """
}
