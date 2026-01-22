process BIN_TAX_REPORT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_single'

    input:
	path(gtdbtk_tax)

    output:
        path("*.png"), emit: gtdb_tax
        path("*.csv"), emit: mag_report

    script:

        """
        Rscript /mnt/Bins_tax.R 
	"""
}
