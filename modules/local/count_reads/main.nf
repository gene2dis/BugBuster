process COUNT_READS {

    label 'process_single'

    input:
        tuple val(meta), path(reads)

    output:
	tuple val(meta), path(reads), emit: reads
        path(reads), emit: reads_coassembly
        path("*_fastp_report.tsv"), emit: reads_report

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        R1_in_count=`zcat ${reads[0]} | wc -l`
        R1_r_count=`echo \$((\${R1_in_count}/4))`
        final_reads_count=`echo \$((\${R1_r_count}*2))`

	if [[ ${reads[2]} != null ]]; then
		Singleton_l_count=`zcat  ${reads[2]} | wc -l`
		Singleton_r_count=`echo \$((\${Singleton_l_count}/4))`
		echo "Id\tRaw reads\tRaw reads singletons" > ${prefix}_fastp_report.tsv
		echo "${prefix}\t\${final_reads_count}\t\${Singleton_r_count}" >> ${prefix}_fastp_report.tsv
	else
		echo "Id\tRaw reads" > ${prefix}_fastp_report.tsv
                echo "${prefix}\t\${final_reads_count}" >> ${prefix}_fastp_report.tsv
	fi
	""" 
}
