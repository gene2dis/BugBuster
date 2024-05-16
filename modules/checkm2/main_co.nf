process CHECKM2 {
    container 'quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0'

    cpus 15

    publishDir "${params.output}/workflow/co_assembly/checkm2", pattern: '*quality_report.tsv'
    publishDir "${params.output}/Reports/bins_quality", mode: 'copy', pattern: '*quality_report.tsv'

    input:
        path(bins)
	path(checkm_db)

    output:
        path("*quality_report.tsv"), emit: all_reports
        path("*_metawrap_quality_report.tsv"), emit: metawrap_report

    script:
        def prefix = "co_assembly"

        """
        checkm2 predict \\
                --threads $task.cpus \\
                --input ${bins[0]} \\
                --output-directory ${prefix}_metabat_checkm2_report \\
                --database_path ${checkm_db} \\
                -x .fa 

        checkm2 predict \\
                --threads $task.cpus \\
                --input ${bins[1]} \\
                --output-directory ${prefix}_semibin_checkm2_report \\
                --database_path ${checkm_db} \\
                -x .fa 

        checkm2 predict \\
                --threads $task.cpus \\
                --input ${bins[2]} \\
                --output-directory ${prefix}_autometa_checkm2_report \\
                --database_path ${checkm_db} \\
                -x .fa 

	 checkm2 predict \\
                --threads $task.cpus \\
                --input ${bins[3]} \\
                --output-directory ${prefix}_metawrap_checkm2_report \\
                --database_path ${checkm_db} \\
                -x .fa

        mv ${prefix}_metabat_checkm2_report/quality_report.tsv ${prefix}_metabat_quality_report.tsv 
        mv ${prefix}_semibin_checkm2_report/quality_report.tsv ${prefix}_semibin_quality_report.tsv 
        mv ${prefix}_autometa_checkm2_report/quality_report.tsv ${prefix}_autometa_quality_report.tsv 
        mv ${prefix}_metawrap_checkm2_report/quality_report.tsv ${prefix}_metawrap_quality_report.tsv
	"""
}
