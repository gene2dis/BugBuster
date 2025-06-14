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
        blobtools create \\
                  --infile ${contigs} \\
                  --hitsfile ${blastn_hits} \\
                  --nodes ${tax_files}/nodes.dmp \\
                  --names ${tax_files}/names.dmp \\
                  --bam ${bam} \\
                  --out ${prefix} 

        blobtools view \\
                  --input ${prefix}.blobDB.json \\
                  --out ${prefix}_Blob_table \\
                  --rank all
	"""
}
