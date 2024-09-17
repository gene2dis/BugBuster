process BIN_QUALITY_REPORT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    cpus 1

    publishDir "${params.output}/Figures/Bins_reports", mode: 'copy', pattern: '*.tiff'
    publishDir "${params.output}/Reports/Bins_quality", mode: 'copy', pattern: '*.csv'

    input:
        path(reports)

    output:
        tuple path("Sample_bins_quality_plot.tiff"), path("Total_bins_quality_plot.tiff"), path("Mag_summary.csv"), path("Mag_table.csv")

    script:
        """
        Rscript /mnt/Bin_checkm_per_sample_plot.R
        Rscript /mnt/Bin_checkm_general_plot.R
	"""
}
