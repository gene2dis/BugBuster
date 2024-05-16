process METABAT2 {
    container 'quay.io/biocontainers/metabat2:2.15--h4da6f23_2'

    cpus 5

    publishDir "${params.output}/workflow/${meta.id}/metabat2", pattern: '*metabat_bins'
    publishDir "${params.output}/Assembly/${meta.id}/Bins", mode: 'copy', pattern: '*metabat_bins'

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
