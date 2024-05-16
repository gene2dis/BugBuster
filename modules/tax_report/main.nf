process TAX_REPORT {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/Reports/Reads_report", mode: 'copy', pattern: 'Reads_report.csv'
    publishDir "${params.output}/Figures/Read_level_taxonomy/Classified_reads", mode: 'copy', pattern: 'Kraken_plot.tiff'

    input:
        path(reports)

    output:
        path("Reads_report.csv"), emit: report
        path("Kraken_plot.tiff"), emit: figure

    script:
        def db_used = "${params.kraken_db_used}"

        """
	Rscript /mnt/Tax_unify_report.R ${db_used}
        mv All_reports.csv Reads_report.csv
        """
}
