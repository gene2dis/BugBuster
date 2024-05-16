process BOWTIE2_SAMTOOLS {
    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    cpus 5

    publishDir "${params.output}/workflow/${meta.id}/bowtie2_samtools", pattern: '*_all_reads.bam'

    input:
        tuple val(meta), path(reads), path(bins)

    // Aqui creo los canales de salida de este proceso.

    output:
        tuple val(meta), path("*_all_reads.bam"), emit: reads

    script:
        def prefix = "${meta.id}"

        """
        for bin in ${bins}/*.fa; do
            bin_name=`echo \${bin} | sed -E 's/.+?b//g' | sed 's/.fa//g' | sed 's/in/bin/g'`

	    bowtie2-build \${bin} ${prefix}_bins_index

            bowtie2 \\
            -x ${prefix}_bins_index \\
            -p $task.cpus \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
	    -S ${prefix}_\${bin_name}_paired_reads.sam 2> ${prefix}_bowtie_map.log 

	    samtools sort -@ $task.cpus \\
                  -o ${prefix}_\${bin_name}_paired_reads.bam \\
                  ${prefix}_\${bin_name}_paired_reads.sam

	    if [[ ${reads[2]} == null ]]; then
		    mv ${prefix}_\${bin_name}_paired_reads.bam ${prefix}_\${bin_name}_all_reads.bam
            fi

	    if [[ ${reads[2]} != null ]]; then
	            bowtie2 \\
	            -x ${prefix}_bins_index \\
	            -p $task.cpus \\
	            -U ${reads[2]} \\
	            -S ${prefix}_\${bin_name}_singletons.sam 2> ${prefix}_bowtie_singleton_map.log

		    samtools sort -@ $task.cpus \\
                          -o ${prefix}_\${bin_name}_singletons_reads.bam \\
                          ${prefix}_\${bin_name}_singletons.sam

		    samtools merge ${prefix}_\${bin_name}_all_reads.bam ${prefix}_\${bin_name}_paired_reads.bam ${prefix}_\${bin_name}_singletons_reads.bam
	            rm -f ${prefix}_bowtie_singleton_map.log
                    rm -f ${prefix}_\${bin_name}_singletons.sam
            fi
            rm -f ${prefix}_\${bin_name}_paired_reads.sam
            rm -f ${prefix}_bins_index*
            rm -f ${prefix}_bowtie_map.log
        done
	""" 
}
