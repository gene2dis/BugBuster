process BOWTIE2 {
    container 'quay.io/biocontainers/bowtie2:2.5.3--py310ha0a81b8_0'

    label 'process_medium'

    publishDir "${params.output}/workflow/${meta.id}/bowtie_${db_alias}", pattern: '*_map_*.fastq.gz'

    input:
        tuple val(meta), path(reads)
	path index_db
	val index_basename
	val db_alias

    output:
        tuple val(meta), path("*map*fastq.gz"), emit: reads
	path("*_bowtie_report.tsv"), emit: report
        path("*map*fastq.gz"), emit: reads_coassembly	

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        bowtie2 \\
        $args \\
        -x ${index_db}/${index_basename} \\
        -p $task.cpus \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --un-conc-gz ${prefix}_${db_alias} > ${prefix}_bowtie_map_${db_alias}

        mv ${prefix}_${db_alias}.1 ${prefix}_R1_map_${db_alias}.fastq.gz &&
        mv ${prefix}_${db_alias}.2 ${prefix}_R2_map_${db_alias}.fastq.gz &&

        R1_in_count=`zcat ${prefix}_R1_map_${db_alias}.fastq.gz | wc -l`
        R2_in_count=`zcat ${prefix}_R2_map_${db_alias}.fastq.gz | wc -l`
        R1_r_count=`echo \$((\${R1_in_count}/4))`
        R2_r_count=`echo \$((\${R2_in_count}/4))`
        final_reads_count=`echo \$((\${R1_r_count}+\${R2_r_count}))`


	if [[ ${reads[2]} != null ]]; then
		bowtie2 \\
                $args \\
		-x ${index_db}/${index_basename} \\
		-p $task.cpus \\
		-U ${reads[2]} \\
		--un-gz ${prefix}_${db_alias} > ${prefix}_bowtie_singleton_map_${db_alias}

		mv ${prefix}_${db_alias} ${prefix}_Singleton_map_${db_alias}.fastq.gz
		rm -f ${prefix}_bowtie_singleton_map_${db_alias}
	
	fi

	if [[ ${reads[2]} != null ]]; then
		Singleton_l_count=`zcat  ${prefix}_Singleton_map_${db_alias}.fastq.gz | wc -l`
		Singleton_r_count=`echo \$((\${Singleton_l_count}/4))`
		echo "Id\tBowtie ${db_alias}\tBowtie ${db_alias} singletons" > ${prefix}_${db_alias}_bowtie_report.tsv
		echo "${prefix}\t\${final_reads_count}\t\${Singleton_r_count}" >> ${prefix}_${db_alias}_bowtie_report.tsv
	else
		echo "Id\tBowtie ${db_alias}" > ${prefix}_${db_alias}_bowtie_report.tsv
                echo "${prefix}\t\${final_reads_count}" >> ${prefix}_${db_alias}_bowtie_report.tsv
	fi
	rm -f ${prefix}_bowtie_map_${db_alias}
	""" 
}
