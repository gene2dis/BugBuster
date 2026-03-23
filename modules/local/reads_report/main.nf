process READS_REPORT {

    container 'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0'

    label 'process_single'

    input:
        path(reports)
        val(args)

    output:
	path("*")

    script:
        """
        report_unify.py ${args}
	""" 
}
