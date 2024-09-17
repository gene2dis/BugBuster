process KRAKEN_TO_PHYLOSEQ {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/Reports/read_taxonomy", mode: 'copy', pattern: '*.RDS'
    publishDir "${params.output}/Figures/Read_level_taxonomy/Tax_plot", mode: 'copy', pattern: '*.tiff'

    cpus 1

    input:
        path(biom)

    output:
        path("*.RDS"), emit: kraken_to_phyloseq
        path("*.tiff"), emit: tiff

    script:

        """
        Rscript /mnt/Tax_kraken_to_phyloseq.R
        """
}
