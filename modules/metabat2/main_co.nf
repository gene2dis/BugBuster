process METABAT2 {
    container 'quay.io/biocontainers/metabat2:2.15--h4da6f23_2'

    cpus 5

    publishDir "${params.output}/workflow/co_assembly/metabat2", pattern: 'metabat_bins'
    publishDir "${params.output}/Co_assembly/Bins", mode: 'copy', pattern: 'metabat_bins'

    input:
        tuple path(contigs), path(depth)

    output:
        path("*metabat_bins"), emit: metabat2

    script:
        """
	metabat2 -i ${contigs} -a ${depth} -o metabat_bins/co_assembly_metabat_bin -t $task.cpus
	"""
}
