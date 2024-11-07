process BOWTIE2_SAMTOOLS {
    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    label 'process_medium'

    publishDir "${params.output}/workflow/${meta.id}/bowtie2_samtools", pattern: '*_all_reads.bam'

    input:
        tuple val(meta), path(reads), path(contigs)

    output:
        tuple val(meta), path(contigs), path("*_all_reads.bam"), emit: contigs_and_bam
        tuple val(meta), path("*_all_reads.bam"), emit: only_bam

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
                rm -f ${prefix}_paired_reads.bam
                rm -f ${prefix}_singletons_reads.bam
	fi
        rm -f ${prefix}_contigs_index*
        rm -f ${prefix}_bowtie_map.log
        rm -f ${prefix}_paired_reads.sam
	""" 
}

process BOWTIE2_SAMTOOLS_COASSEMBLY {
    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    label 'process_medium'

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
                rm -f ${prefix}_paired_reads.bam
                rm -f ${prefix}_singletons_reads.bam
                rm -f ${prefix}_singletons.sam
                rm -f ${prefix}_bowtie_singleton_map.log
	fi
        rm -f ${prefix}_bowtie_map.log
        rm -f ${prefix}_paired_reads.sam
        rm -f co_assembly_contigs_index*
	""" 
}

process BOWTIE2_SAMTOOLS_DEPTH {
    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    label 'process_medium'

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
