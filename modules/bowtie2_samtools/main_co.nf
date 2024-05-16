process BOWTIE2_SAMTOOLS {
    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    cpus 20

    publishDir "${params.output}/workflow/co_assembly/bowtie2_samtools", pattern: 'co_assembly_all_reads.bam'

    input:
        tuple val(meta), path(reads), path(contigs)

    // Aqui creo los canales de salida de este proceso.

    output:
        path("*_all_reads.bam"), emit: bam

    script:
        def prefix = "${meta.id}"

        """
        R1_list=`ls | grep -E '.+?R1_map.+' | tr '\\n' ',' | sed 's/.\$//'`
        R2_list=`ls | grep -E '.+?R2_map.+' | tr '\\n' ',' | sed 's/.\$//'`
        Singleton_list=`ls | grep -E '.+?Singleton_map.+' | tr '\\n' ',' | sed 's/.\$//'`

	bowtie2-build ${contigs} co_assembly_contigs_index

        bowtie2 \\
        -x co_assembly_contigs_index \\
        -p $task.cpus \\
        -1 \$R1_list \\
        -2 \$R2_list \\
	-S ${prefix}_paired_reads.sam 2> ${prefix}_bowtie_map.log 

	samtools sort -@ $task.cpus \\
              -o ${prefix}_paired_reads.bam \\
              ${prefix}_paired_reads.sam

        if [[ `ls | grep 'Singleton_map'` == "" ]]; then
		mv ${prefix}_paired_reads.bam ${prefix}_all_reads.bam
        fi

        if [[ `ls | grep 'Singleton_map'` != "" ]]; then
	        bowtie2 \\
	        -x co_assembly_contigs_index \\
	        -p $task.cpus \\
	        -U \$Singleton_list \\
	        -S ${prefix}_singletons.sam 2> ${prefix}_bowtie_singleton_map.log

		samtools sort -@ $task.cpus \\
                      -o ${prefix}_singletons_reads.bam \\
                      ${prefix}_singletons.sam

		samtools merge ${prefix}_all_reads.bam ${prefix}_paired_reads.bam ${prefix}_singletons_reads.bam
                rm -f ${prefix}_singletons.sam
                rm -f ${prefix}_bowtie_singleton_map.log
	fi
        rm -f ${prefix}_bowtie_map.log
        rm -f ${prefix}_paired_reads.sam
        rm -f co_assembly_contigs_index*
	""" 
}
