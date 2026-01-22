process BIN_QUALITY_REPORT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_single'

    input:
        path(reports)

    output:
        tuple path("*.png"), path("*.csv")

    script:
        """
        Rscript /mnt/Bin_checkm_general_plot.R
	"""
}
