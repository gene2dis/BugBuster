process SAMTOOLS_INDEX {
    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    cpus 1

    publishDir "${params.output}/workflow/${meta.id}/samtool_index", pattern: '*_all_reads.bam'

    input:
        tuple val(meta), path(sorted_bam)

    // Aqui creo los canales de salida de este proceso.

    output:
        tuple val(meta), path(sorted_bam), path("*.bam.bai")

    script:
        def prefix = "${meta.id}"

        """
	samtools index ${sorted_bam} -@ $task.cpus
	""" 
}
