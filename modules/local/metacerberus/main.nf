process METACERBERUS_CONTIGS {
    container 'quay.io/ffuentessantander/metacerberus:1.2.1'

    label 'process_high'

    publishDir "${params.output}/workflow/${meta.id}/metacerberus_contigs", pattern: '*_annotation'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/contigs", mode: 'copy', pattern: '*_annotation_results'

    input:
        tuple val(meta), path(contigs)

    output:
        tuple val(meta), path("*annotation"), emit: reads
        path("*_annotation_results"), emit: results

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        metacerberus.py \\
                     --prodigal ${contigs} \\
                     ${args} \\
                     --meta \\
                     --cpus $task.cpus \\
                     --dir_out ${prefix}_annotation

        mv ${prefix}_annotation/step_10-visualizeData ${prefix}_annotation_results
        """
}

process METACERBERUS_BINS {
    container 'quay.io/ffuentessantander/metacerberus:1.2.1'

    label 'process_high'

    publishDir "${params.output}/workflow/${meta.id}/metacerberus_bins", pattern: '*_annotation'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/bins", mode: 'copy', pattern: '*_annotation_results'

    input:
        tuple val(meta), path(bin_folder)

    output:
        tuple val(meta), path("*annotation"), emit: reads
        path("*_annotation_results"), emit: results

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        for bin in ${bin_folder}/*; do
            bin_name=`echo \$bin | sed 's/.fa//g' | sed 's/./_/g'`
            metacerberus.py \\
                        --prodigal \${bin} \\
                        ${args} \\
                        --cpus $task.cpus \\
                        --dir_out ${prefix}_\${bin_name}_annotation

            mv ${prefix}_\${bin_name}_annotation/step_10-visualizeData ${prefix}_\${bin_name}_annotation_results
        done
        """
}
