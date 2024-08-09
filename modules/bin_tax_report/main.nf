process BIN_TAX_REPORT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    cpus 1

    publishDir "${params.output}/Figures/Bins_reports", mode: 'copy', pattern: 'Bins_tax_plot.tiff'

    input:
	path(gtdbtk_tax)

    output:
        path("Bins_tax_plot.tiff"), emit: gtdb_tax

    script:

        """
        Rscript /mnt/Bins_tax.R
	"""
}
