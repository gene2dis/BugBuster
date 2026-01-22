process ARG_BLOBPLOT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_single'

    input:
        path(blob_table)

    output:
        path("*.png")

    script:
        """
        Rscript /mnt/ARG_blob_plot.R
	"""
}
