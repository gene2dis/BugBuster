process METABULI_CLASSIFY {
    tag "${meta.id}"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/metabuli:1.2.0--pl5321h0bb26bb_0' :
        'quay.io/biocontainers/metabuli:1.2.0--pl5321h0bb26bb_0' }"

    label 'process_medium'

    input:
        tuple val(meta), path(contigs)
        path(db)

    output:
        tuple val(meta), path("*_classifications.tsv"), emit: classifications
        tuple val(meta), path("*_metabuli_hits.tsv")  , emit: hits
        path("*_report.tsv")                          , emit: report
        path "versions.yml"                           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        // --max-ram caps Metabuli's own buffers; 90 % of the task's memory
        // leaves room for the rest of the process (32 GiB of process_medium's
        // 36 GB, the Phase 9c setting: 5-12 GB used, same calls at 32-128)
        def max_ram = Math.max(1, (task.memory.toGiga() * 0.9) as int)
        """
        # --seq-mode 3 (long reads) is Metabuli's single-sequence mode for
        # contigs; --lineage 1 adds the full lineage that
        # ARG_CONTIG_LEVEL_REPORT reads
        metabuli classify \\
            --seq-mode 3 \\
            --lineage 1 \\
            --threads ${task.cpus} \\
            --max-ram ${max_ram} \\
            ${args} \\
            ${contigs} \\
            ${db} \\
            metabuli_out \\
            ${prefix}

        mv metabuli_out/${prefix}_classifications.tsv metabuli_out/${prefix}_report.tsv .

        # BlobTools hits file (contig, taxid, score): classified contigs only,
        # the others become BlobTools' no-hit
        awk -F '\\t' '!/^#/ && \$1 == "1" { print \$2 "\\t" \$3 "\\t" \$5 }' \\
            ${prefix}_classifications.tsv > ${prefix}_metabuli_hits.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            metabuli: \$(metabuli 2>&1 | sed -n 's/^metabuli Version: //p')
            metabuli_db: created \$(sed -n 's/^Creation_date\\t//p' ${db}/db.parameters)
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"
        """
        printf '#is_classified\\tname\\ttaxID\\tquery_length\\tscore\\te_value\\trank\\tlineage\\ttaxID:match_count\\n' > ${prefix}_classifications.tsv
        touch ${prefix}_metabuli_hits.tsv ${prefix}_report.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            metabuli: 1.2.0
            metabuli_db: created 2026-4-3
        END_VERSIONS
        """
}
