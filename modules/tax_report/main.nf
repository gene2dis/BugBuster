process TAX_REPORT {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/reports/read_level/", mode: 'copy', pattern: 'reads_summary.csv'
    publishDir "${params.output}/reports/read_level/taxonomy", mode: 'copy', pattern: 'classified_reads.tiff'

    label 'process_single'

    input:
        path(reports)

    output:
        path("reads_summary.csv"), emit: report
        path("classified_reads.tiff"), emit: figure

    script:
        def db_used = "${params.kraken_db_used}"

        """
	Rscript /mnt/Tax_unify_report.R ${db_used}
        mv All_reports.csv reads_summary.csv
        mv Kraken_plot.tiff classified_reads.tiff
        """
}
