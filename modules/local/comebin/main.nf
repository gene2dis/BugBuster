process COMEBIN {
    container 'quay.io/biocontainers/comebin:1.0.4--hdfd78af_0'

    label 'process_high'

    publishDir "${params.output}/workflow/${meta.id}/comebin", pattern: '*_comebin_bins/comebin_res/comebin_res_bins'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/raw_bins", mode: 'copy', pattern: '*_comebin_bins/comebin_res/comebin_res_bins'

    input:
        tuple val(meta), path(contigs), path(bam)
    
    output:
        tuple val(meta), path("*_comebin_bins/comebin_res/comebin_res_bins"), emit: comebin

    script:
        def prefix = "${meta.id}"

        """
        mkdir -p ${prefix}_comebin_bins/comebin_res/comebin_res_bins

        contig_count=\$(grep -c "^>" ${contigs} || true)
        
        if [[ \$contig_count -lt 10 ]]; then
            echo "Warning: Sample ${prefix} has only \$contig_count contigs. Skipping COMEBIN (requires at least 10 contigs)." >&2
            echo "COMEBIN skipped due to insufficient contigs (\$contig_count < 10)" > ${prefix}_comebin_bins/comebin_res/comebin_res_bins/SKIPPED.txt
            exit 0
        fi

        set +e
        if [[ $task.attempt == 1 ]]; then
               bash run_comebin.sh \\
                   -a ${contigs} \\
                   -p ./ \\
                   -o ${prefix}_comebin_bins \\
                   -n 8 \\
                   -t $task.cpus
               exit_code=\$?
        elif [[ $task.attempt == 2 ]]; then 
                bash run_comebin.sh \\
                   -a ${contigs} \\
                   -p ./ \\
                   -o ${prefix}_comebin_bins \\
                   -n 8 \\
                   -t $task.cpus \\
                   -b 896 \\
                   -e 1792 \\
                   -c 1792
                exit_code=\$?
        elif [[ $task.attempt == 3 ]]; then
                bash run_comebin.sh \\
                   -a ${contigs} \\
                   -p ./ \\
                   -o ${prefix}_comebin_bins \\
                   -n 8 \\
                   -t $task.cpus \\
                   -b 512 \\
                   -e 1024 \\
                   -c 1024
                exit_code=\$?
        fi
        set -e

        if [[ \$exit_code -ne 0 ]]; then
            echo "Warning: COMEBIN failed for sample ${prefix}. Creating empty bin directory to allow pipeline to continue." >&2
            mkdir -p ${prefix}_comebin_bins/comebin_res/comebin_res_bins
            echo "COMEBIN failed - insufficient data or features for binning" > ${prefix}_comebin_bins/comebin_res/comebin_res_bins/FAILED.txt
            exit 0
        fi

	"""
}

process COMEBIN_COASSEMBLY {
    container 'quay.io/biocontainers/comebin:1.0.4--hdfd78af_0'

    label 'process_high'

    publishDir "${params.output}/workflow/co_assembly/comebin", pattern: 'coassembly_comebin_bins/comebin_res/comebin_res_bins'
    publishDir "${params.output}/contigs_and_bins/co_assembly/raw_bins", mode: 'copy', pattern: 'coassembly_comebin_bins/comebin_res/comebin_res_bins'

    input:
        path(bams_and_contigs)

    output:
        path("coassembly_comebin_bins/comebin_res/comebin_res_bins"), emit: comebin

    script:
        def prefix = "coassembly"

        """
        contigs=`ls | grep -E '.+?filtered_contigs.+' | tr '\\n' ',' | sed 's/.\$//'`

        mkdir -p ${prefix}_comebin_bins/comebin_res/comebin_res_bins

        contig_count=0
        IFS=',' read -ra CONTIG_FILES <<< "\${contigs}"
        for contig_file in "\${CONTIG_FILES[@]}"; do
            count=\$(grep -c "^>" "\$contig_file" || true)
            contig_count=\$((contig_count + count))
        done
        
        if [[ \$contig_count -lt 10 ]]; then
            echo "Warning: Coassembly has only \$contig_count contigs. Skipping COMEBIN (requires at least 10 contigs)." >&2
            echo "COMEBIN skipped due to insufficient contigs (\$contig_count < 10)" > ${prefix}_comebin_bins/comebin_res/comebin_res_bins/SKIPPED.txt
            exit 0
        fi

        set +e
        if [[ $task.attempt == 1 ]]; then
               bash run_comebin.sh \\
                   -a \${contigs} \\
                   -p ./ \\
                   -o ${prefix}_comebin_bins \\
                   -n 8 \\
                   -t $task.cpus
               exit_code=\$?
        elif [[ $task.attempt == 2 ]]; then
                bash run_comebin.sh \\
                   -a \${contigs} \\
                   -p ./ \\
                   -o ${prefix}_comebin_bins \\
                   -n 8 \\
                   -t $task.cpus \\
                   -b 896 \\
                   -e 1792 \\
                   -c 1792
                exit_code=\$?
        elif [[ $task.attempt == 3 ]]; then
                bash run_comebin.sh \\
                   -a \${contigs} \\
                   -p ./ \\
                   -o ${prefix}_comebin_bins \\
                   -n 8 \\
                   -t $task.cpus \\
                   -b 512 \\
                   -e 1024 \\
                   -c 1024
                exit_code=\$?
        fi
        set -e

        if [[ \$exit_code -ne 0 ]]; then
            echo "Warning: COMEBIN failed for coassembly. Creating empty bin directory to allow pipeline to continue." >&2
            mkdir -p ${prefix}_comebin_bins/comebin_res/comebin_res_bins
            echo "COMEBIN failed - insufficient data or features for binning" > ${prefix}_comebin_bins/comebin_res/comebin_res_bins/FAILED.txt
            exit 0
        fi

        """
}
