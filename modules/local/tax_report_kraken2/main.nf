process TAX_REPORT_KRAKEN2 {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/reports/read_level/", mode: 'copy', pattern: '*.csv', saveAs: { filename -> "Reads_report.csv" }
    publishDir "${params.output}/reports/read_level/taxonomy", mode: 'copy', pattern: '*.png'

    label 'process_single'

    input:
        path(reports)

    output:
        path("*.csv"), emit: report
        path(reports), emit: input
        path("*.png"), emit: figure

    script:
        def db_used = "${params.kraken_db_used}"

        """
	Rscript /mnt/Tax_unify_report.R ${db_used}
        cp -L Reads_report.csv tmp_report.csv
        """
}
