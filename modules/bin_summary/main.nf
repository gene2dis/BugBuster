process BIN_SUMMARY {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/reports", mode: 'copy', pattern: 'Bin_summary.csv'

    label 'process_single'

    input:
        path(reports)

    output:
        path("Bin_summary.csv")

    script:
        """
	Rscript /mnt/Bin_summary.R
	"""
}
