process BLOBPLOT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    cpus 1

    publishDir "${params.output}/Figures", mode: 'copy', pattern: 'Blob_plot.tiff'

    input:
        path(blob_table)

    output:
        path("Blob_plot.tiff")

    script:
        """
        Rscript /mnt/Blobplot.R 
	"""
}
