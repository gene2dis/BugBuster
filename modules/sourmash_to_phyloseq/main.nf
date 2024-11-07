process SOURMASH_TO_PHYLOSEQ {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_single'

    publishDir "${params.output}/reports/read_level/taxonomy", mode: 'copy', pattern: '*.RDS'
    publishDir "${params.output}/reports/read_level/taxonomy", mode: 'copy', pattern: '*.png'

    input:
        path(sourmash_gather)
	val db_name

    output:
        path("*.RDS"), emit: ps_object
	path("*.png"), emit: plots

    script:

        """
        Rscript /mnt/Tax_sourmash_to_phyloseq.R ${db_name}
        """
}
