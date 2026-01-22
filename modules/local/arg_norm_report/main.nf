process ARG_NORM_REPORT {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_single'

    input:
        path(arg_data)

    output:
        path("*.csv"), emit: arg_reports

    script: 

    """
    Rscript /mnt/Read_arg_norm.R
    """

}
