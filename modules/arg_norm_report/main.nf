process ARG_NORM_REPORT {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    publishDir "${params.output}/reports/read_level/arg_prediction", mode: 'copy', pattern: '*.csv'

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
