process BLOBTOOLS {
    container 'quay.io/ffuentessantander/blobtools:1.1.1'

    tag "${meta.id}"
    label 'process_single'

    input:
        tuple val(meta), path(contigs), path(blastn_hits), path(bam), path(bam_bai), path(tax_files)

    output:
        tuple val(meta), path("*_Blob_tabl*"), emit: blob_table
        path("*_Blob_tabl*"), emit: only_blob
        path "versions.yml", emit: versions

    script:
        def prefix = "${meta.id}"

        """
        # tax_files stages either loose taxdump files (FORMAT_TAXDUMP_FILES)
        # or a directory (--custom_taxdump_files) — locate them either way
        nodes_dmp=\$(find -L . -maxdepth 2 -name nodes.dmp | head -n 1)
        names_dmp=\$(find -L . -maxdepth 2 -name names.dmp | head -n 1)

        # --db: write the parsed nodesDB into the task dir; the default is
        # inside the blobtools package, which is read-only in a
        # singularity/apptainer image (design doc Q21)
        blobtools create \\
                  --infile ${contigs} \\
                  --hitsfile ${blastn_hits} \\
                  --nodes \${nodes_dmp} \\
                  --names \${names_dmp} \\
                  --db nodesDB.txt \\
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
