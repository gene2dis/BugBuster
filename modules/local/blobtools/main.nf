process BLOBTOOLS {
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/blobtools:1.1.1--py_1' :
        'quay.io/biocontainers/blobtools:1.1.1--py_1' }"

    tag "${meta.id}"
    label 'process_single'

    input:
        tuple val(meta), path(contigs), path(hits), path(bam), path(bam_bai)
        path(nodesdb)

    output:
        tuple val(meta), path("*_Blob_tabl*"), emit: blob_table
        path("*_Blob_tabl*"), emit: only_blob
        path "versions.yml", emit: versions

    script:
        def prefix = "${meta.id}"

        """
        # nodesdb: BLOBTOOLS_NODESDB's parse of the Metabuli database's
        # taxdump, built once for every sample; read here, not written
        blobtools create \\
                  --infile ${contigs} \\
                  --hitsfile ${hits} \\
                  --db ${nodesdb} \\
                  --bam ${bam} \\
                  --out ${prefix}

        blobtools view \\
                  --input ${prefix}.blobDB.json \\
                  --out ${prefix}_Blob_table \\
                  --rank all

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blobtools: \$(blobtools --version 2>&1 | tail -1 || echo 1.1.1)
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"

        """
        touch ${prefix}_Blob_table.blobDB.table.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blobtools: 1.1.1
        END_VERSIONS
        """
}
