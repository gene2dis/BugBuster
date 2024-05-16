process KARGA {
    container 'quay.io/ffuentessantander/karga:1.0'

    cpus 1

    publishDir "${params.output}/workflow/${meta.id}/karga", pattern: '*KARGA_mappedGenes.csv'

    input:
        tuple val(meta), path(report), path(reads)
	path(params.karga_db)	    

    output:
        tuple val(meta), path("*_all_reads_KARGA_mappedGenes.csv"), emit: kargva

    script:
        def prefix = "${meta.id}"

        """
	mv /mnt/* .
        java KARGA k:17 d:${params.karga_db} r:n -Xmx32GB ${reads[0]}
        """
}
