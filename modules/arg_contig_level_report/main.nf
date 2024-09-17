process ARG_CONTIG_LEVEL_REPORT {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/ARG_prediction/Contig_level", mode: 'copy', pattern: 'Contig_tax_and_arg_prediction.tsv'

    cpus 1

    input:
        path(arg_contig_data)

    output:
        path("Contig_tax_and_arg_prediction.tsv"), emit: arg_reports

    script: 

    """
    Rscript /mnt/Contig_arg_unify.R  
    """

}
