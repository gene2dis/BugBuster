process KRAKEN_TO_PHYLOSEQ {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/reports/read_level/taxonomy", mode: 'copy', pattern: '*.RDS'

    label 'process_single'

    input:
        path(biom)

    output:
        path("*.RDS"), emit: kraken_to_phyloseq

    script:

        """
        Rscript /mnt/Tax_kraken_to_phyloseq.R
        """
}
