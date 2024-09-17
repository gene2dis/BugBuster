process BLOBTOOLS {
    container 'quay.io/biocontainers/blobtools:1.1.1--py_1'

    cpus 1

    input:
        tuple val(meta), path(contigs), path(blastn_hits), path(bam), path(bam_bai)
        path(NCBI_NODES_DMP)
        path(NCBI_NAMES_DMP)

    output:
        tuple val(meta), path("*_Blob_tabl*"), emit: blob_table
        path("*_Blob_tabl*"), emit: only_blob

    script:
        def prefix = "${meta.id}"

        """
        blobtools create \\
                  --infile ${contigs} \\
                  --hitsfile ${blastn_hits} \\
                  --nodes ${NCBI_NODES_DMP} \\
                  --names ${NCBI_NAMES_DMP} \\
                  --bam ${bam} \\
                  --out ${prefix} 

        blobtools view \\
                  --input ${prefix}.blobDB.json \\
                  --out ${prefix}_Blob_table \\
                  --rank all
	"""
}
