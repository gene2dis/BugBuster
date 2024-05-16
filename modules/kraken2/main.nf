process KRAKEN2 {
    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0'

    maxForks 4
    cpus 5

    publishDir "${params.output}/workflow/${meta.id}/kraken2_${db_name}", mode: 'copy', pattern: '*_kraken_report'

    input:
        tuple val(meta), path(reads)
	path db
	val db_name

    output:
        tuple val(meta), path("*_kraken_report"), emit: kraken
	path("*_report.tsv"), emit: report

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        kraken2 \\
                --db $db \\
                --threads $task.cpus \\
                --report ${prefix}_${db_name}_kraken_report \\
                --gzip-compressed \\
                $args \\
                --paired ${reads[0]} ${reads[1]} \\
                --output ${prefix}_${db_name}_kraken_output 
	
	echo "Id\tKraken DB\tUnclassified\tClassified" > ${meta.id}_${db_name}_report.tsv
        head -1 ${prefix}_${db_name}_kraken_report | awk -v id=${prefix} -v db=${db_name} '{print id "\t" db "\t" \$1 "\t" 100 - \$1}' >> ${meta.id}_${db_name}_report.tsv

        """
}
