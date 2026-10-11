process NT_BLASTN {
    container 'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'

    // Full-nt megablast: memory-hungry (label scales memory/time by attempt
    // under the global retry-on-OOM errorStrategy)
    tag "${meta.id}"
    label 'process_high'

    input:
        tuple val(meta), path(contigs), path(nt_db)

    output:
        tuple val(meta), path(contigs), path("*_megablast.out"), emit: megablast_to_blob
        path "versions.yml", emit: versions

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        # The DB is staged either as a directory (--custom_blast_db) or as
        # loose volume files (FORMAT_NT_BLAST_DB). Resolve the DB prefix from
        # the alias file (multi-volume) or a single .nin index.
        nal_file=\$(find -L . -maxdepth 2 -name '*.nal' | head -n 1)
        if [ -n "\${nal_file}" ]; then
            db_prefix="\${nal_file%.nal}"
        else
            nin_file=\$(find -L . -maxdepth 2 -name '*.nin' | head -n 1)
            db_prefix=\$(echo "\${nin_file}" | sed -E 's/(\\.[0-9]+)?\\.nin\$//')
        fi
        if [ -z "\${db_prefix}" ]; then
            echo "ERROR: no BLAST database (*.nal/*.nin) found among staged inputs" >&2
            exit 1
        fi

        blastn \\
            -task megablast \\
            -query ${contigs} \\
            -db \${db_prefix} \\
            -outfmt '6 qseqid staxids bitscore std' \\
            -num_threads $task.cpus \\
            ${args} \\
            -out ${prefix}_assembly_vs_nt_megablast.out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blastn: \$(blastn -version 2>&1 | head -n 1 | sed 's/^.*blastn: //')
        END_VERSIONS
	"""

    stub:
        def prefix = "${meta.id}"

        """
        touch ${prefix}_assembly_vs_nt_megablast.out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blastn: 2.15.0
        END_VERSIONS
        """
}
