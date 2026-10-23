process DEEPARG_BINS {
    container 'quay.io/ffuentessantander/deeparg:1.0.4'
        
    tag "${meta.id}"
    label 'process_medium'

    input:
        tuple val(meta), path(prodigal_bins), path(deeparg_db)

    output:
        tuple val(meta), path("*_deeparg_results"), emit: deeparg_bins
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        # The proteins dir may be empty (upstream bins dir held only a
        # SKIPPED/FAILED marker) — nullglob keeps the *.faa loop from running
        # deeparg on a literal '*.faa'
        shopt -s nullglob

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

        # The image's cp is BusyBox: without -L it copies the staged input
        # directory as a symlink, and ../ inside it is no longer the task dir
        cp -rL ${prodigal_bins} tmp_bins
        mkdir ${prefix}_deeparg_results
        cd tmp_bins

        # Up to ${task.cpus} bins run at once; their threads share the task's CPUs
        bins=(*.faa)
        jobs=\$(( \${#bins[@]} < ${task.cpus} ? \${#bins[@]} : ${task.cpus} ))
        export DEEPARG_DIAMOND_THREADS=\$(( ${task.cpus} / (jobs > 0 ? jobs : 1) ))
        export OPENBLAS_NUM_THREADS=\$DEEPARG_DIAMOND_THREADS OMP_NUM_THREADS=\$DEEPARG_DIAMOND_THREADS

        # Parallel processing of bins using background jobs
        pids=()
        for bin_prot in *.faa; do
            bin_name=\$(echo \$bin_prot | sed -E 's/_proteins.faa//g')
            
            # Run deeparg in background
            (
                deeparg predict \\
                    -d ../${deeparg_db} \\
                    --model LS \\
                    --type prot \\
                    $args \\
                    --input \${bin_prot} \\
                    --out \${bin_name}_deep_arg.out
            ) &
            pids+=(\$!)
            
            # When we reach max concurrent jobs, wait for each to finish
            if [ \${#pids[@]} -ge ${task.cpus} ]; then
                for pid in "\${pids[@]}"; do wait "\$pid" || exit 1; done
                pids=()
            fi
        done
        
        # Wait for all remaining jobs
        for pid in "\${pids[@]}"; do wait "\$pid" || exit 1; done
        
        for arg_file in *.mapping.ARG; do mv "\$arg_file" ../${prefix}_deeparg_results/; done
        cd ..
        rm -rf tmp_bins

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeparg: \$(pip show deeparg 2>/dev/null | sed -n 's/Version: //p' || echo 1.0.4)
        END_VERSIONS
	"""

    stub:
        def prefix = "${meta.id}"

        """
        mkdir ${prefix}_deeparg_results

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeparg: 1.0.4
        END_VERSIONS
        """
}

process DEEPARG_CONTIGS {
    container 'quay.io/ffuentessantander/deeparg:1.0.4'

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

        deeparg predict \\
                -d ${deeparg_db} \\
                --model LS \\
                --type prot \\
                $args \\
                --input ${prefix}_proteins_input.faa \\
                --out ${prefix}_contigs_deep_arg.out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeparg: \$(pip show deeparg 2>/dev/null | sed -n 's/Version: //p' || echo 1.0.4)
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
