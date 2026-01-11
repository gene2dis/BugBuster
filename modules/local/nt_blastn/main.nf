process NT_BLASTN {
    container 'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'

    label 'process_low'

    input:
        tuple val(meta), path(contigs), path(nt_db)

    output:
        tuple val(meta), path(contigs), path("*_megablast.out"), emit: megablast_to_blob

    script:
        def prefix = "${meta.id}"

        """
	db_basename=`ls ${nt_db}/ | grep -E '\\.[0]+?\\.' | cut -f1 -d'.' | sort | uniq`

        blastn \\
            -task megablast \\
            -query ${contigs} \\
            -db ${nt_db}/\${db_basename} \\
            -outfmt '6 qseqid staxids bitscore std' \\
            -max_target_seqs 1 \\
            -max_hsps 1 \\
            -num_threads $task.cpus \\
            -evalue 1e-25 \\
            -out ${prefix}_assembly_vs_nt_megablast.out
	"""
}
