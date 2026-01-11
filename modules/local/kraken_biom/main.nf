process KRAKEN_BIOM {
    container 'quay.io/biocontainers/kraken-biom:1.2.0--pyh5e36f6f_0'

    maxForks 4 
    label 'process_single'

    publishDir "${params.output}/reports/read_level/taxonomy", mode: 'copy', pattern: '*.biom'

    input:
        path(bracken_report)
        val db_name

    output:
        path("*.biom"), emit: bracken_report

    script:

        """
	kraken-biom *.report --fmt json -o ${db_name}_microbiom.biom
        """
}
