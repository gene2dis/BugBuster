process RGI_BWT {
    container 'quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0'
    
    label 'process_medium'
    
    tag "${meta.id}"

    input:
        tuple val(meta), path(reads)
        path(card_db)

    output:
        tuple val(meta), path("*.allele_mapping_data.txt"), emit: allele_mapping
        tuple val(meta), path("*.gene_mapping_data.txt"), emit: gene_mapping
        tuple val(meta), path("*.sorted.length_100.bam"), emit: bam
        path("*.overall_mapping_stats.txt"), emit: overall_stats
        path("*.reference_mapping_stats.txt"), emit: reference_stats
        path("*.artifacts_mapping_stats.txt"), emit: artifacts_stats, optional: true
        path("versions.yml"), emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        def aligner = params.rgi_aligner ?: 'kma'
        def wildcard = params.rgi_include_wildcard ? '--include_wildcard' : ''
        
        // Handle paired-end reads with optional singletons
        def read_one = reads[0]
        def read_two = reads.size() > 1 ? reads[1] : ''
        
        """
        # Copy database to local directory for RGI access
        cp -r ${card_db} rgi_db
        export RGI_DATA=\$(pwd)/rgi_db
        
        # Run RGI bwt
        rgi bwt \\
            --read_one ${read_one} \\
            ${read_two ? "--read_two ${read_two}" : ""} \\
            --aligner ${aligner} \\
            --threads ${task.cpus} \\
            --output_file ${prefix}_rgi_bwt \\
            --local \\
            ${wildcard} \\
            ${args}

        # Verify output files were created
        if [ ! -f "${prefix}_rgi_bwt.allele_mapping_data.txt" ]; then
            echo "ERROR: RGI bwt did not produce expected output files"
            exit 1
        fi

        # Create versions file
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: \$(rgi --version 2>&1 | grep -oP 'rgi \\K[0-9.]+' || echo "unknown")
            ${aligner}: \$(${aligner} --version 2>&1 | head -1 || echo "unknown")
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"
        """
        touch ${prefix}_rgi_bwt.allele_mapping_data.txt
        touch ${prefix}_rgi_bwt.gene_mapping_data.txt
        touch ${prefix}_rgi_bwt.sorted.length_100.bam
        touch ${prefix}_rgi_bwt.overall_mapping_stats.txt
        touch ${prefix}_rgi_bwt.reference_mapping_stats.txt
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: 6.0.3
            kma: 1.4.14
        END_VERSIONS
        """
}
