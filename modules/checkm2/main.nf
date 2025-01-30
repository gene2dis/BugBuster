process CHECKM2_COASSEMBLY {
    container 'quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0'

    label 'process_medium'

    publishDir "${params.output}/workflow/co_assembly/checkm2", pattern: '*quality_report.tsv'

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

process CHECKM2 {
    container 'quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0'

    label 'process_medium'

    publishDir "${params.output}/workflow/${meta.id}/checkm2", pattern: '*quality_report.tsv'

    input:
        tuple val(meta), path(metabat2), path(semibin), path(comebin), path(metawrap)
	path(checkm_db)

    output:
        path("*quality_report.tsv"), emit: all_reports
        path("*_metawrap_quality_report.tsv"), emit: metawrap_report

    script:
        def prefix = "${meta.id}"

        """
        checkm2 predict \\
                --threads $task.cpus \\
                --input ${metabat2} \\
                --output-directory ${prefix}_metabat_checkm2_report \\
                --database_path ${checkm_db} \\
                -x .fa 

        checkm2 predict \\
                --threads $task.cpus \\
                --input ${semibin} \\
                --output-directory ${prefix}_semibin_checkm2_report \\
                --database_path ${checkm_db} \\
                -x .fa 

        checkm2 predict \\
                --threads $task.cpus \\
                --input ${comebin} \\
                --output-directory ${prefix}_comebin_checkm2_report \\
                --database_path ${checkm_db} \\
                -x .fa 

	 checkm2 predict \\
                --threads $task.cpus \\
                --input ${metawrap} \\
                --output-directory ${prefix}_metawrap_checkm2_report \\
                --database_path ${checkm_db} \\
                -x .fa

        mv ${prefix}_metabat_checkm2_report/quality_report.tsv ${prefix}_metabat_quality_report.tsv 
        mv ${prefix}_semibin_checkm2_report/quality_report.tsv ${prefix}_semibin_quality_report.tsv 
        mv ${prefix}_comebin_checkm2_report/quality_report.tsv ${prefix}_comebin_quality_report.tsv 
        mv ${prefix}_metawrap_checkm2_report/quality_report.tsv ${prefix}_metawrap_quality_report.tsv
	"""
}
