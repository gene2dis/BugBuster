process BLOBOLOGY_COASSEMBLY {
    container 'quay.io/biocontainers/metawrap:1.2--hdfd78af_2'

    cpus 40

    publishDir "${params.output}/workflow/co_assembly/blobology", pattern: '*metawrap*bins'
    publishDir "${params.output}/Co_assembly/", mode: 'copy', pattern: '*blobology'

    input:
        path(contigs)
        path(bins)
	path(reads)

    output:
        path("*blobology"), emit: blobology

    script:
        def prefix = "co_assembly"
        """
        metawrap blobology \\
                 -a ${contigs} \\
                 -t $task.cpus \\
                 -o ${prefix}_blobology \\
                 --bins ${bins} \\
                 ${reads}_fastq
	"""
}
