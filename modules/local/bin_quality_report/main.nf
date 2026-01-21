process BIN_QUALITY_REPORT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_single'

    publishDir "${params.output}/reports/bin_level/bins_quality", mode: 'copy', pattern: '*.png'
    publishDir "${params.output}/reports/bin_level/bins_quality", mode: 'copy', pattern: '*.csv'

    input:
        path(reports)

    output:
        tuple path("*.png"), path("*.csv")

    script:
        """
        Rscript /mnt/Bin_checkm_general_plot.R
	"""
}
