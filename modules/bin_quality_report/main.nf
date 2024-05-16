process BIN_QUALITY_REPORT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    cpus 15

    publishDir "${params.output}/Figures/Bins_reports", mode: 'copy', pattern: 'Bins_quality_plot.tiff'

    input:
        path(reports)

    output:
        path("Bins_quality_plot.tiff")

    script:
        """
	Rscript /mnt/Checkm_unify.R 
	"""
}
