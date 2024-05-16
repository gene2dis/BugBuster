process QFILTER {

    publishDir "${params.output}/Filtered_reads/qreport", mode: 'copy', pattern: '*.txt'

    input:
        tuple val(meta), path(reads), path(json)

    // Aqui creo los canales de salida de este proceso.

    output:
        tuple val(meta), path(reads), emit: qfilter

    shell:
        def prefix = "${meta.id}"

        """

	""" 
}
