process GTDB_TK {
    container 'quay.io/ffuentessantander/gtdbtk:2.4.0'

    label 'process_high'
    maxForks 1

    publishDir "${params.output}/workflow/${meta.id}/gtdb-tk", pattern: '*_gtdbtk_*'

    input:
        tuple val(meta), path(metawrap)
	path(gtdbtk_db)

    output:
        tuple val(meta), path("*_gtdbtk_*"), emit: gtdb_tk
        path("*_gtdbtk_*"), emit: report

    script:
        def prefix = "${meta.id}"

        """
        export GTDBTK_DATA_PATH="${gtdbtk_db}"

        gtdbtk classify_wf \\
               --genome_dir ${metawrap} \\
               --out_dir ${prefix}_gtdb \\
               --mash_db ${gtdbtk_db}/mash \\
               --cpus $task.cpus \\
               --extension .fa \\
               --pplacer_cpus $task.cpus
        
	mv ${prefix}_gtdb/classify/gtdbtk.bac120.summary.tsv ${prefix}_gtdbtk_bac120.tsv

	if [ -f ${prefix}_gtdb/classify/gtdbtk.ar53.summary.tsv ]; then 
                mv ${prefix}_gtdb/classify/gtdbtk.ar53.summary.tsv ${prefix}_gtdbtk_ar53.tsv
	fi

	"""
}

process GTDB_TK_COASSEMBLY {
    container 'quay.io/ffuentessantander/gtdbtk:2.4.0'

    label 'process_high'

    publishDir "${params.output}/workflow/co_assembly/gtdb-tk", pattern: '*_gtdbtk_*'

    input:
        path(metawrap)
	path(gtdbtk_db)

    output:
        path("*_gtdbtk_*"), emit: gtdb_tk
        path("*_gtdbtk_*"), emit: report        

    script:
        def prefix = "co_assembly"

        """
        export GTDBTK_DATA_PATH="${gtdbtk_db}"

        gtdbtk classify_wf \\
               --genome_dir ${metawrap} \\
               --out_dir ${prefix}_gtdb \\
               --mash_db ${gtdbtk_db}/mash \\
               --cpus $task.cpus \\
               --extension .fa \\
               --pplacer_cpus $task.cpus
        
	mv ${prefix}_gtdb/classify/gtdbtk.bac120.summary.tsv ${prefix}_gtdbtk_bac120.tsv

	if [ -f ${prefix}_gtdb/classify/gtdbtk.ar53.summary.tsv ]; then 
                mv ${prefix}_gtdb/classify/gtdbtk.ar53.summary.tsv ${prefix}_gtdbtk_ar53.tsv
	fi

	"""
}
