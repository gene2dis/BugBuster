process RGI_KMER {
    container 'quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0'
    
    label 'process_medium'
    
    tag "${meta.id}"

    input:
        tuple val(meta), path(bam)
        path(card_db)

    output:
        tuple val(meta), path("*_61mer_analysis.json"), emit: kmer_json
        tuple val(meta), path("*_61mer_analysis.txt"), emit: kmer_txt
        tuple val(meta), path("*.allele.txt"), emit: allele_predictions, optional: true
        tuple val(meta), path("*.gene.txt"), emit: gene_predictions, optional: true
        path("versions.yml"), emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        def kmer_size = params.rgi_kmer_size ?: 61
        def minimum = params.rgi_min_kmer_coverage ?: 10
        
        """
        # Copy database to local directory for RGI access
        cp -r ${card_db} rgi_db
        export RGI_DATA=\$(pwd)/rgi_db
        
        # Run RGI kmer_query for pathogen-of-origin prediction
        rgi kmer_query \\
            --bwt \\
            --rgi_bwt_bam ${bam} \\
            --kmer_size ${kmer_size} \\
            --threads ${task.cpus} \\
            --minimum ${minimum} \\
            --output ${prefix}_rgi_kmer \\
            --local \\
            ${args}

        # Verify output files were created
        if [ ! -f "${prefix}_rgi_kmer_${kmer_size}mer_analysis.json" ]; then
            echo "WARNING: RGI kmer_query did not produce expected JSON output"
            # Create empty output to prevent pipeline failure
            echo '{}' > ${prefix}_rgi_kmer_${kmer_size}mer_analysis.json
            echo "No k-mer results" > ${prefix}_rgi_kmer_${kmer_size}mer_analysis.txt
        fi

        # Create versions file
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: \$(rgi --version 2>&1 | grep -oP 'rgi \\K[0-9.]+' || echo "unknown")
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"
        def kmer_size = params.rgi_kmer_size ?: 61
        """
        touch ${prefix}_rgi_kmer_${kmer_size}mer_analysis.json
        touch ${prefix}_rgi_kmer_${kmer_size}mer_analysis.txt
        touch ${prefix}_rgi_kmer.allele.txt
        touch ${prefix}_rgi_kmer.gene.txt
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: 6.0.3
        END_VERSIONS
        """
}
