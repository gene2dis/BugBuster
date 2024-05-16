process BLOBOLOGY {
    container 'quay.io/biocontainers/metawrap:1.2--hdfd78af_2'

    cpus 10

    publishDir "${params.output}/workflow/${meta.id}/blobology", pattern: '*_blobology'
    publishDir "${params.output}/Assembly/${meta.id}/Refined_bins/blobology", mode: 'copy', pattern: '*_blobology'

    input:
        tuple val(meta), path(contigs), path(bins), path(reads)

    output:
        tuple val(meta), path("*_blobology"), emit: metawrap

    script:
        def prefix = "${meta.id}"

        """
        metawrap blobology \\
                 -a ${contigs} \\
                 -t $task.cpus \\
                 -o ${prefix}_blobology \\
                 --bins ${bins} \\
                 ${reads}_fastq
	"""
}
