process DEEPARG_BINS {
    container 'quay.io/ffuentessantander/deeparg:1.0.4'
        
    label 'process_medium'

    input:
        tuple val(meta), path(prodigal_bins), path(deeparg_db)

    output:
        tuple val(meta), path("*_deeparg_results"), emit: deeparg_bins

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        cp -r ${prodigal_bins} tmp_bins
        mkdir ${prefix}_deeparg_results
        cd tmp_bins
        
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
        
        mv *.mapping.ARG ../${prefix}_deeparg_results/
        cd ..
        rm -rf tmp_bins
	"""
}

process DEEPARG_CONTIGS {
    container 'quay.io/ffuentessantander/deeparg:1.0.4'

    label 'process_medium'

    input:
        tuple val(meta), path(prodigal_contigs), path(deeparg_db)

    output:
        tuple val(meta), path("*.ARG"), emit: deeparg
        path("*.ARG"), emit: only_deeparg

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        deeparg predict \\
                -d ${deeparg_db} \\
                --model LS \\
                --type prot \\
                $args \\
                --input ${prodigal_contigs} \\
                --out ${prefix}_contigs_deep_arg.out
        """
}
