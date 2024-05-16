process MEGAHIT {
    container 'quay.io/biocontainers/megahit:1.2.9--h43eeafb_4'

    cpus 5

    publishDir "${params.output}/workflow/${meta.id}/megahit", pattern: '*_contigs.fa'

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
