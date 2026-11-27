process DEEPARG_CONTIGS {
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/deeparg:1.0.4--pyhdfd78af_0' :
        'quay.io/biocontainers/deeparg:1.0.4--pyhdfd78af_0' }"

    tag "${meta.id}"
    label 'process_medium'

    input:
        tuple val(meta), path(prodigal_contigs), path(deeparg_db)

    output:
        tuple val(meta), path("*.ARG"), emit: deeparg
        path("*.ARG"), emit: only_deeparg
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        # nf-core prodigal emits gzipped proteins; deeparg cannot read gzip.
        # gzip -cdf also passes plain FASTA through unchanged.
        gzip -cdf ${prodigal_contigs} > ${prefix}_proteins_input.faa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeparg: \$(v=\$(pip show deeparg 2>/dev/null | sed -n 's/Version: //p'); echo \${v:-1.0.4})
        END_VERSIONS

        # An assembly whose contigs give Pyrodigal no gene (design doc Q24, option
        # b; roadmap Phase 29, CL5): deeparg predict stops on an empty FASTA
        # (IOError ... align.daa.tsv). Write the two tables a run without hits
        # writes, DeepARG 1.0.4's header row alone, without running it
        if ! grep -q '^>' ${prefix}_proteins_input.faa; then
            echo "WARNING: ${prefix} has no predicted protein: DeepARG skipped, header-only ARG tables written" >&2
            for table in mapping.ARG mapping.potential.ARG; do
                printf '#ARG\\tquery-start\\tquery-end\\tread_id\\tpredicted_ARG-class\\tbest-hit\\tprobability\\tidentity\\talignment-length\\talignment-bitscore\\talignment-evalue\\tcounts\\n' \
                    > ${prefix}_contigs_deep_arg.out.\${table}
            done
            rm ${prefix}_proteins_input.faa
            exit 0
        fi

        # deeparg predict runs `diamond blastp` with no thread count, so DIAMOND
        # starts one thread per host core whatever the task's CPUs, and the
        # image's OpenBLAS does the same in the classifier. A diamond wrapper
        # first on PATH adds --threads. deeparg ignores DIAMOND's exit status
        # (a DIAMOND killed for memory gave an empty ARG table and a "successful"
        # task), so the wrapper records a failure and the task exits with
        # DIAMOND's status after deeparg (137 is retried with more memory)
        DEEPARG_REAL_DIAMOND=\$(command -v diamond) || { echo "diamond not found" >&2; exit 1; }
        export DEEPARG_REAL_DIAMOND
        export DEEPARG_DIAMOND_FAILED="\$PWD/diamond_failed"
        mkdir diamond_wrapper
        cat > diamond_wrapper/diamond <<-'END_WRAPPER'
        #!/bin/sh
        case "\$1" in
            blastp|blastx|view)
                "\$DEEPARG_REAL_DIAMOND" "\$@" --threads "\$DEEPARG_DIAMOND_THREADS"
                status=\$?
                # keep the first failure: the later 'view' then fails too
                [ \$status -eq 0 ] || [ -s "\$DEEPARG_DIAMOND_FAILED" ] || echo \$status > "\$DEEPARG_DIAMOND_FAILED"
                exit \$status ;;
        esac
        exec "\$DEEPARG_REAL_DIAMOND" "\$@"
        END_WRAPPER
        chmod +x diamond_wrapper/diamond
        export PATH="\$PWD/diamond_wrapper:\$PATH"
        export DEEPARG_DIAMOND_THREADS=${task.cpus}
        export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus}
        # Theano creates its compile cache under \$HOME/.theano: '/' in a
        # docker task without HOME (not writable), the user's home under
        # singularity. Keep it in the task directory (design doc Q21)
        export THEANO_FLAGS="base_compiledir=\$PWD/theano_compiledir"

        deeparg_status=0
        deeparg predict \\
                -d ${deeparg_db} \\
                --model LS \\
                --type prot \\
                $args \\
                --input ${prefix}_proteins_input.faa \\
                --out ${prefix}_contigs_deep_arg.out || deeparg_status=\$?

        # DIAMOND's status first: deeparg then fails (or not) for lack of its
        # alignments, and only DIAMOND's 137 says "retry with more memory"
        if [ -s diamond_failed ]; then
            echo "ERROR: diamond exited with status \$(cat diamond_failed) inside deeparg predict; its ARG table would be incomplete" >&2
            exit \$(cat diamond_failed)
        fi
        [ \$deeparg_status -eq 0 ] || exit \$deeparg_status
        """

    stub:
        def prefix = "${meta.id}"

        """
        touch ${prefix}_contigs_deep_arg.out.mapping.ARG

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeparg: 1.0.4
        END_VERSIONS
        """
}
