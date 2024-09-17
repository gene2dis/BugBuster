process BIN_TAX_REPORT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    cpus 1

    publishDir "${params.output}/Figures/Bins_reports", mode: 'copy', pattern: 'MAGs_tax_plot.tiff'
    publishDir "${params.output}/Reports/Bins_taxonomy", mode: 'copy', pattern: 'MAGs_tax_summary.csv'

    input:
	path(gtdbtk_tax)

    output:
        path("MAGs_tax_plot.tiff"), emit: gtdb_tax
        path("MAGs_tax_summary.csv"), emit: mag_report

    script:

        """
        cp /mnt/Bins_tax.R .
        Rscript /mnt/Bins_tax.R 
	"""
}
