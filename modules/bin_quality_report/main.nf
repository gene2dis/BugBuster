process BIN_QUALITY_REPORT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_single'

    publishDir "${params.output}/reports/bins_quality", mode: 'copy', pattern: '*.tiff'
    publishDir "${params.output}/reports/bins_quality", mode: 'copy', pattern: '*.csv'

    input:
        path(reports)

    output:
        tuple path("Total_bins_quality_plot.tiff"), path("Mag_quality_summary.csv"), path("All_Mag_quality_table.csv"), path("Refined_Mag_quality_table.csv")

    script:
        """
        Rscript /mnt/Bin_checkm_general_plot.R
	"""
}
