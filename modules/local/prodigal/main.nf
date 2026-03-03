process PRODIGAL_BINS {
    container 'quay.io/biocontainers/prodigal:2.6.3--h031d066_8'

    label 'process_medium'

    input:
        tuple val(meta), path(metawrap)

    output:
        tuple val(meta), path("*_bins_proteins"), emit: prodigal_bins

    script:
        def prefix = "${meta.id}"

        """
        #!/bin/bash
        set -euo pipefail

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
        
        mv *_genes.gff ${prefix}_bins_genes
        mv *_proteins.faa ${prefix}_bins_proteins
        cd ..
        mv tmp_bins/${prefix}_bins_genes/ .
        mv tmp_bins/${prefix}_bins_proteins/ .

        rm -rf refined_bins/
	"""
}

process PRODIGAL_CONTIGS {
    container 'quay.io/biocontainers/prodigal:2.6.3--h031d066_8'

    label 'process_single'

    input:
        tuple val(meta), path(contigs)

    output:
        tuple val(meta), path("*_contigs_proteins.faa"), emit: prodigal_contigs

    script:
        def prefix = "${meta.id}"

        """
        prodigal -i ${contigs} \\
                 -o ${prefix}_contigs_genes.gff \\
                 -a ${prefix}_contigs_proteins.faa \\
                 -p meta 
        """
}
