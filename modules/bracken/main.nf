process BRACKEN {
    container 'quay.io/biocontainers/bracken:2.9--py39h1f90b4d_0'

    maxForks 4
    label 'process_single'

    input:
        tuple val(meta), path(kraken_report)
        path db
        val db_name

    output:
        path("*.report"), emit: bracken_report

    script:
        def prefix = "${meta.id}"
        def read_len = "${params.bracken_read_len}"
        def tax_level = "${params.bracken_tax_level}"

        """
        bracken \\
                -d ${db} \\
                -i ${kraken_report} \\
                -o ${prefix}_${db_name}.braken \\
                -w ${prefix}_${db_name}.report \\
                -r $read_len \\
                -l $tax_level
        """
}
