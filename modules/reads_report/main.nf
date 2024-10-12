process READS_REPORT {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/reports/read_level", mode: 'copy', pattern: 'Box_plot_reads.tiff'

    label 'process_single'

    input:
        path(reports)

    output:
	path("All_reads_report.csv"), emit: report
	path("Box_plot_reads.tiff"), emit: figure

    script:
        def db_used = "${params.bowtie_db_used}"

        """
        Rscript /mnt/Report_unify.R $db_used 
	""" 
}
