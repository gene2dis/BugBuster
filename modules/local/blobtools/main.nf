process BLOBTOOLS {
    container 'quay.io/ffuentessantander/blobtools:1.1.1'

    label 'process_single'

    input:
        tuple val(meta), path(contigs), path(blastn_hits), path(bam), path(bam_bai), path(tax_files)

    output:
        tuple val(meta), path("*_Blob_tabl*"), emit: blob_table
        path("*_Blob_tabl*"), emit: only_blob

    script:
        def prefix = "${meta.id}"

        """
        # tax_files stages either loose taxdump files (FORMAT_TAXDUMP_FILES)
        # or a directory (--custom_taxdump_files) — locate them either way
        nodes_dmp=\$(find -L . -maxdepth 2 -name nodes.dmp | head -n 1)
        names_dmp=\$(find -L . -maxdepth 2 -name names.dmp | head -n 1)

        blobtools create \\
                  --infile ${contigs} \\
                  --hitsfile ${blastn_hits} \\
                  --nodes \${nodes_dmp} \\
                  --names \${names_dmp} \\
                  --bam ${bam} \\
                  --out ${prefix}

        blobtools view \\
                  --input ${prefix}.blobDB.json \\
                  --out ${prefix}_Blob_table \\
                  --rank all
	"""

    stub:
        def prefix = "${meta.id}"

        """
        touch ${prefix}_Blob_table.blobDB.table.txt
        """
}
