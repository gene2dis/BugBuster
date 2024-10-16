process BLOBPLOT {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_single'

    publishDir "${params.output}/reports/contig_level", mode: 'copy', pattern: 'Phylum_blob_plot.tiff'

    input:
        path(blob_table)

    output:
        path("Phylum_blob_plot.tiff")

    script:
        """
        Rscript /mnt/Blobplot.R 
	"""
}