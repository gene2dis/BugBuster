process ARG_CONTIG_LEVEL_REPORT {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_single'

    input:
        path(arg_contig_data)

    output:
        path("Contig_tax_and_arg_prediction.tsv"), emit: arg_reports

    script: 

    """
    Rscript /mnt/Contig_arg_unify.R  
    """

}
