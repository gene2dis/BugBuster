process COUNT_READS {

    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'ubuntu:22.04' }"

    input:
        tuple val(meta), path(reads)

    output:
	tuple val(meta), path(reads), emit: reads
        path(reads), emit: reads_coassembly
        path("*_fastp_report.tsv"), emit: reads_report
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def _args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        R1_in_count=`zcat -f ${reads[0]} | wc -l`
        R1_r_count=`echo \$((\${R1_in_count}/4))`
        final_reads_count=`echo \$((\${R1_r_count}*2))`

	if [[ ${reads[2]} != null ]]; then
		Singleton_l_count=`zcat -f ${reads[2]} | wc -l`
		Singleton_r_count=`echo \$((\${Singleton_l_count}/4))`
		printf 'Id\\tRaw reads\\tRaw reads singletons\\n' > ${prefix}_fastp_report.tsv
		printf '%s\\t%s\\t%s\\n' "${prefix}" "\${final_reads_count}" "\${Singleton_r_count}" >> ${prefix}_fastp_report.tsv
	else
		printf 'Id\\tRaw reads\\n' > ${prefix}_fastp_report.tsv
		printf '%s\\t%s\\n' "${prefix}" "\${final_reads_count}" >> ${prefix}_fastp_report.tsv
	fi

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    gzip: \$(gzip --version | head -n 1 | awk '{print \$NF}')
	END_VERSIONS
	"""

    stub:
        def prefix = "${meta.id}"

        """
        echo -e "Id\\tRaw reads" > ${prefix}_fastp_report.tsv
        echo -e "${prefix}\\t1000" >> ${prefix}_fastp_report.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gzip: 1.10
        END_VERSIONS
        """
}
