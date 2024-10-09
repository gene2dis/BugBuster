process METABAT2 {
    container 'quay.io/biocontainers/metabat2:2.15--h4da6f23_2'

    label 'process_medium'

    publishDir "${params.output}/workflow/${meta.id}/metabat2", pattern: '*metabat_bins'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/raw_bins", mode: 'copy', pattern: '*metabat_bins'

    input:
        tuple val(meta), path(contigs), path(depth)

    output:
        tuple val(meta), path("*metabat_bins"), emit: metabat2

    script:
        def prefix = "${meta.id}"

        """
	metabat2 -i ${contigs} -a ${depth} -o ${prefix}_metabat_bins/${prefix}_metabat_bin -t $task.cpus
	"""
}

process METABAT2_COASSEMBLY {
    container 'quay.io/biocontainers/metabat2:2.15--h4da6f23_2'

    label 'process_medium'

    publishDir "${params.output}/workflow/co_assembly/metabat2", pattern: 'metabat_bins'
    publishDir "${params.output}/coassembly_and_bins/raw_bins", mode: 'copy', pattern: 'metabat_bins'

    input:
        tuple path(contigs), path(depth)

    output:
        path("*metabat_bins"), emit: metabat2

    script:
        """
	metabat2 -i ${contigs} -a ${depth} -o metabat_bins/co_assembly_metabat_bin -t $task.cpus
	"""
}
