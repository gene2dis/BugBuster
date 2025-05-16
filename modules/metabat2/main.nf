process METABAT2 {
    container 'quay.io/biocontainers/metabat2:2.15--h4da6f23_2'

    label 'process_medium'

    publishDir "${params.output}/workflow/${meta.id}/metabat2", pattern: '*metabat_bins'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/raw_bins", mode: 'copy', pattern: '*metabat_bins'

    input:
        tuple val(meta), path(contigs), path(depth)

    output:
        tuple val(meta), path("*metabat_bins"), emit: metabat2

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
	metabat2 -i ${contigs} \\
		 -a ${depth} \\
		 -o ${prefix}_metabat_bins/${prefix}_metabat_bin \\
		 $args \\
		 -t $task.cpus
	"""
}

process METABAT2_COASSEMBLY {
    container 'quay.io/biocontainers/metabat2:2.15--h4da6f23_2'

    label 'process_medium'

    publishDir "${params.output}/workflow/co_assembly/metabat2", pattern: 'metabat_bins'
    publishDir "${params.output}/contigs_and_bins/co_assembly/raw_bins", mode: 'copy', pattern: 'metabat_bins'

    input:
        tuple path(contigs), path(depth)

    output:
        path("*metabat_bins"), emit: metabat2

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''

        """
	metabat2 -i ${contigs} \\
		 -a ${depth} \\
		 -o metabat_bins/co_assembly_metabat_bin \\
                 $args \\
		 -t $task.cpus
	"""
}
