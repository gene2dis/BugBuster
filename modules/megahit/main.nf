process MEGAHIT_COASSEMBLY {
    container 'quay.io/biocontainers/megahit:1.2.9--h43eeafb_4'

    label 'process_high'

    input:
        path(reads)
    
    output:
        tuple path(reads), path("co_assembly_contigs.fa"), emit: megahit

    script:

        """
        R1_list=`ls | grep -E '.+?R1.+' | tr '\\n' ',' | sed 's/.\$//'`
        R2_list=`ls | grep -E '.+?R2.+' | tr '\\n' ',' | sed 's/.\$//'`
        Singleton_list=`ls | grep -E '.+?Singleton.+' | tr '\\n' ',' | sed 's/.\$//'`

	if [[ `ls | grep 'Singleton_map'` != "" ]]; then	
		megahit \\
                   -1 \$R1_list \\
                   -2 \$R2_list \\
                   --read \$Singleton_list \\
                   -t $task.cpus \\
                   -o co_assembly
	else
		megahit \\
                -1 \$R1_list \\
                -2 \$R2_list \\
                -t $task.cpus \\
                -o co_assembly
	fi
	
	mv co_assembly/final.contigs.fa co_assembly_contigs.fa
	
	"""
}

process MEGAHIT {
    container 'quay.io/biocontainers/megahit:1.2.9--h43eeafb_4'

    label 'process_medium'

    input:
        tuple val(meta), path(reads)
    
    output:
        tuple val(meta), path(reads), path("*_contigs.fa"), emit: contigs_and_reads

    script:
        def prefix = "${meta.id}"

        """
	if [[ ${reads[2]} != null ]]; then	
		megahit \\
                   -1 ${reads[0]} \\
                   -2 ${reads[1]} \\
                   --read ${reads[2]} \\
                   -t $task.cpus \\
                   -o ${prefix}_assembly
	else
		megahit \\
                -1 ${reads[0]} \\
                -2 ${reads[1]} \\
                -t $task.cpus \\
                -o ${prefix}_assembly
	fi
	
	mv ${prefix}_assembly/final.contigs.fa ${prefix}_contigs.fa
	
	"""
}
