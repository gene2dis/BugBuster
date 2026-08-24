process METACERBERUS_CONTIGS {
    container 'quay.io/ffuentessantander/metacerberus:1.2.1'

    label 'process_high'

    input:
        tuple val(meta), path(contigs)

    output:
        tuple val(meta), path("*annotation"), emit: reads
        path("*_annotation_results"), emit: results
        path "versions.yml", emit: versions

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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            metacerberus: \$(metacerberus.py --version 2>&1 | grep -o '[0-9.]*' | head -n 1 || echo 1.2.1)
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"

        """
        mkdir ${prefix}_annotation
        mkdir ${prefix}_annotation_results

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            metacerberus: 1.2.1
        END_VERSIONS
        """
}

process METACERBERUS_BINS {
    container 'quay.io/ffuentessantander/metacerberus:1.2.1'

    label 'process_high'

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
