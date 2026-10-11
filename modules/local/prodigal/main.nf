process PRODIGAL_BINS {
    container 'quay.io/biocontainers/prodigal:2.6.3--h031d066_8'

    tag "${meta.id}"
    label 'process_medium'

    input:
        tuple val(meta), path(metawrap)

    output:
        tuple val(meta), path("*_bins_proteins"), emit: prodigal_bins
        path "versions.yml", emit: versions

    script:
        def prefix = "${meta.id}"

        """
        #!/bin/bash
        set -euo pipefail
        # A bins dir may hold only a SKIPPED/FAILED marker (see METAWRAP) —
        # without nullglob the *.fa loop would feed prodigal a literal '*.fa'
        shopt -s nullglob

        cp -r ${metawrap} tmp_bins
        cd tmp_bins

        mkdir ${prefix}_bins_genes
        mkdir ${prefix}_bins_proteins

        # Improved parallelization: run up to task.cpus jobs concurrently
        pids=()
        for file in *.fa; do
            file_name=\$(echo \$file | sed 's/.fa//g')
            
            # Run prodigal in background
            (
                prodigal -i \$file \\
                    -o ${prefix}_\${file_name}_genes.gff \\
                    -a ${prefix}_\${file_name}_proteins.faa \\
                    -p single
            ) &
            pids+=(\$!)
            
            # When we reach max concurrent jobs, wait for each to finish
            if [ \${#pids[@]} -ge ${task.cpus} ]; then
                for pid in "\${pids[@]}"; do wait "\$pid" || exit 1; done
                pids=()
            fi
        done
        
        # Wait for all remaining jobs to complete
        for pid in "\${pids[@]}"; do wait "\$pid" || exit 1; done
        
        for gff in *_genes.gff; do mv "\$gff" ${prefix}_bins_genes/; done
        for faa in *_proteins.faa; do mv "\$faa" ${prefix}_bins_proteins/; done
        cd ..
        mv tmp_bins/${prefix}_bins_genes/ .
        mv tmp_bins/${prefix}_bins_proteins/ .

        rm -rf refined_bins/

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            prodigal: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
        END_VERSIONS
	"""

    stub:
        def prefix = "${meta.id}"

        """
        mkdir ${prefix}_bins_proteins

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            prodigal: 2.6.3
        END_VERSIONS
        """
}
