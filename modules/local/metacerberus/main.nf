process METACERBERUS_READS {
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metacerberus:1.4.0--pyhdfd78af_1' :
        'quay.io/biocontainers/metacerberus:1.4.0--pyhdfd78af_1' }"

    label 'process_high'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.id}_annotation"), emit: annotation
        path("${meta.id}_annotation_results"),          emit: results
        path("versions.yml"),                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        metacerberus.py \\
            --illumina ${reads[0]} ${reads[1]} \\
            --skip-decon \\
            --meta \\
            --cpus $task.cpus \\
            --dir_out ${prefix}_annotation \\
            ${args}

        mv ${prefix}_annotation/step_10-visualizeData ${prefix}_annotation_results

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            metacerberus: \$(metacerberus.py --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        END_VERSIONS
        """
}

process METACERBERUS_CONTIGS {
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metacerberus:1.4.0--pyhdfd78af_1' :
        'quay.io/biocontainers/metacerberus:1.4.0--pyhdfd78af_1' }"

    label 'process_high'

    input:
        tuple val(meta), path(contigs)

    output:
        tuple val(meta), path("${meta.id}_annotation"), emit: annotation
        path("${meta.id}_annotation_results"),          emit: results
        path("versions.yml"),                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args   = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        metacerberus.py \\
            --prodigal ${contigs} \\
            --meta \\
            --cpus $task.cpus \\
            --dir_out ${prefix}_annotation \\
            ${args}

        mv ${prefix}_annotation/step_10-visualizeData ${prefix}_annotation_results

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            metacerberus: \$(metacerberus.py --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        END_VERSIONS
        """
}

process METACERBERUS_BINS_BATCH {
    tag "batch_all_samples"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metacerberus:1.4.0--pyhdfd78af_1' :
        'quay.io/biocontainers/metacerberus:1.4.0--pyhdfd78af_1' }"

    label 'process_high'

    input:
        tuple val(meta_list), path(all_bin_dirs, stageAs: 'bin_?/*')

    output:
        tuple val(meta_list), path("*_annotation_results"), emit: results
        path("versions.yml"),                               emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args     = task.ext.args ?: ''
        def meta_ids = meta_list.collect { m -> m instanceof Map ? (m.id ?: m['id']) : m.toString() }.unique().join(' ')

        """
        set -euo pipefail

        for staged_dir in bin_*/; do
            staged_dir=\${staged_dir%/}
            bin_dir=\$(find "\$staged_dir" -mindepth 1 -maxdepth 1 \\( -type d -o -type l \\) | head -1)
            [ -z "\$bin_dir" ] && continue

            bin_dir_name=\$(basename "\$bin_dir")

            # Extract sample_id from directory name (pattern: {sample_id}_*_bins)
            if [[ "\$bin_dir_name" =~ ^(.+)_(metabat|semibin|comebin|metawrap).*_bins\$ ]]; then
                sample_id="\${BASH_REMATCH[1]}"
            else
                echo "Skipping unrecognized bin directory: \$bin_dir_name"
                continue
            fi

            echo "Processing bins for sample: \$sample_id from \$bin_dir"

            find -L "\$bin_dir" -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) | while read bin_file; do
                bin_name=\$(basename "\$bin_file" | sed 's/\\.[^.]*\$//')
                out_dir="\${sample_id}_\${bin_name}_annotation"

                metacerberus.py \\
                    --prodigal "\$bin_file" \\
                    --meta \\
                    --cpus ${task.cpus} \\
                    --dir_out "\$out_dir" \\
                    ${args}

                if [ -d "\${out_dir}/step_10-visualizeData" ]; then
                    mv "\${out_dir}/step_10-visualizeData" "\${sample_id}_\${bin_name}_annotation_results"
                fi
            done
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            metacerberus: \$(metacerberus.py --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        END_VERSIONS
        """
}
