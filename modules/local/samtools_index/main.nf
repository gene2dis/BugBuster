process SAMTOOLS_INDEX {
    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    cpus 1

    input:
        tuple val(meta), path(sorted_bam)

    output:
        tuple val(meta), path(sorted_bam), path("*.bam.bai")

    script:
        def prefix = "${meta.id}"

        """
	samtools index ${sorted_bam} -@ $task.cpus
	""" 
}
