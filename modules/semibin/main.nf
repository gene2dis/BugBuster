process SEMIBIN {
    container 'quay.io/biocontainers/semibin:2.1.0--pyhdfd78af_0'

    cpus 5

    publishDir "${params.output}/workflow/${meta.id}/semibin", pattern: '*semibin_bins'
    publishDir "${params.output}/Assembly/${meta.id}/Bins", mode: 'copy', pattern: '*semibin_bins'

    input:
        tuple val(meta), path(contigs), path(bam)
    
    output:
        tuple val(meta), path("*semibin_bins/output_bins"), emit: rosella

    script:
        def prefix = "${meta.id}"
        def env_model = "${params.semibin_env_model}"

        """
	SemiBin2 single_easy_bin \\
                 -i ${contigs} \\
                 -b ${bam} \\
                 -o ${prefix}_semibin_bins \\
                 --environment ${env_model} \\
                 --compression none \\
                 --threads $task.cpus
	"""
}
