process READS_REPORT {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/reports/read_level", mode: 'copy', pattern: '*.png'

    label 'process_single'

    input:
        path(reports)

    output:
	path("All_reads_report.csv"), emit: report
	path("*.png"), emit: figure

    script:
        def db_used = "${params.bowtie_db_used}"

        """
        Rscript /mnt/Report_unify.R $db_used 
	""" 
}
