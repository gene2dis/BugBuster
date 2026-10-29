process BIN_ARG_SPLIT {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    tuple val(meta), path(contig_bins), path(proteins), path(deeparg)

    output:
    tuple val(meta), path("*_bins_proteins"), path("*_deeparg_results"), emit: bin_args
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # The bins' ARGs are DEEPARG_CONTIGS' hits on the whole assembly, split by
    # bin (performance review M8): DeepARG min-max scales over its input
    # file, so per-bin runs gave probabilities that depended on the bin
    bin_arg_split.py ${meta.id} ${contig_bins} ${proteins} ${deeparg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${meta.id}_bins_proteins ${meta.id}_deeparg_results
    printf '>k141_1_1\\nMK\\n' > ${meta.id}_bins_proteins/${meta.id}_bin1_proteins.faa
    printf '#ARG\\tquery-start\\tquery-end\\tread_id\\tpredicted_ARG-class\\tbest-hit\\tprobability\\tidentity\\talignment-length\\talignment-bitscore\\talignment-evalue\\tcounts\\n' \\
        > ${meta.id}_deeparg_results/${meta.id}_bin1_deep_arg.out.mapping.ARG
    printf 'ACRA\\t1\\t2\\tk141_1_1\\tmultidrug\\tX|FEATURES|acrA|multidrug|acrA\\t0.99\\t100.0\\t2\\t4.0\\t1e-10\\t1\\n' \\
        >> ${meta.id}_deeparg_results/${meta.id}_bin1_deep_arg.out.mapping.ARG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.12
    END_VERSIONS
    """
}
