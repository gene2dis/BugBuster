process READS_REPORT {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/reports/read_level", mode: 'copy', pattern: '*'

    label 'process_single'

    input:
        path(reports)
        val(args)

    output:
	path("*")

    script:
        """
        Rscript /mnt/Report_unify.R ${args}
	""" 
}
