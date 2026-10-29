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

        # deeparg predict runs `diamond blastp` with no thread count, so DIAMOND
        # starts one thread per host core whatever the task's CPUs, and the
        # image's OpenBLAS does the same in the classifier. A diamond wrapper
        # first on PATH adds --threads; deeparg ignores DIAMOND's exit status,
        # so the real binary must be found here.
        DEEPARG_REAL_DIAMOND=\$(command -v diamond) || { echo "diamond not found" >&2; exit 1; }
        export DEEPARG_REAL_DIAMOND
        mkdir diamond_wrapper
        cat > diamond_wrapper/diamond <<-'END_WRAPPER'
        #!/bin/sh
        case "\$1" in
            blastp|blastx|view) exec "\$DEEPARG_REAL_DIAMOND" "\$@" --threads "\$DEEPARG_DIAMOND_THREADS" ;;
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

        deeparg predict \\
                -d ${deeparg_db} \\
                --model LS \\
                --type prot \\
                $args \\
                --input ${prefix}_proteins_input.faa \\
                --out ${prefix}_contigs_deep_arg.out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeparg: \$(v=\$(pip show deeparg 2>/dev/null | sed -n 's/Version: //p'); echo \${v:-1.0.4})
        END_VERSIONS
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
