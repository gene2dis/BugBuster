process METACERBERUS_CONTIGS {
    container 'quay.io/ffuentessantander/metacerberus:1.2.1'

    label 'process_high'

    publishDir "${params.output}/workflow/${meta.id}/metacerberus_contigs", pattern: '*_annotation'
    publishDir "${params.output}/Annotation/${meta.id}/metacerberus_contigs", mode: 'copy', pattern: '*_annotation_results'

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
