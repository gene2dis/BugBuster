process BIN_TAX_REPORT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_single'

    publishDir "${params.output}/reports/bins_taxonomy", mode: 'copy', pattern: 'MAGs_tax_plot.tiff'
    publishDir "${params.output}/reports/bins_taxonomy", mode: 'copy', pattern: 'MAGs_tax_summary.csv'

    input:
	path(gtdbtk_tax)

    output:
        path("MAGs_tax_plot.tiff"), emit: gtdb_tax
        path("MAGs_tax_summary.csv"), emit: mag_report

    script:

        """
        Rscript /mnt/Bins_tax.R 
	"""
}
