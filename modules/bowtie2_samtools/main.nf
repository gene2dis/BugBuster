process BOWTIE2_SAMTOOLS {
    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    cpus 10

    publishDir "${params.output}/workflow/${meta.id}/bowtie2_samtools", pattern: '*_all_reads.bam'

    input:
        tuple val(meta), path(reads), path(contigs)

    // Aqui creo los canales de salida de este proceso.

    output:
        tuple val(meta), path(contigs), path("*_all_reads.bam"), emit: reads

    script:
        def prefix = "${meta.id}"

        """
	bowtie2-build ${contigs} ${prefix}_contigs_index

        bowtie2 \\
        -x ${prefix}_contigs_index \\
        -p $task.cpus \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
	-S ${prefix}_paired_reads.sam 2> ${prefix}_bowtie_map.log 

	samtools sort -@ $task.cpus \\
              -o ${prefix}_paired_reads.bam \\
              ${prefix}_paired_reads.sam

	if [[ ${reads[2]} == null ]]; then
		mv ${prefix}_paired_reads.bam ${prefix}_all_reads.bam
        fi

	if [[ ${reads[2]} != null ]]; then
	        bowtie2 \\
	        -x ${prefix}_contigs_index \\
	        -p $task.cpus \\
	        -U ${reads[2]} \\
	        -S ${prefix}_singletons.sam 2> ${prefix}_bowtie_singleton_map.log

		samtools sort -@ $task.cpus \\
                      -o ${prefix}_singletons_reads.bam \\
                      ${prefix}_singletons.sam

		samtools merge ${prefix}_all_reads.bam ${prefix}_paired_reads.bam ${prefix}_singletons_reads.bam
                rm -f ${prefix}_singletons.sam
                rm -f ${prefix}_bowtie_singleton_map.log
	fi
        rm -f ${prefix}_contigs_index*
        rm -f ${prefix}_bowtie_map.log
        rm -f ${prefix}_paired_reads.sam
	""" 
}
