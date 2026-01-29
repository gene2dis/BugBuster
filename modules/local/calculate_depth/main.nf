/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALCULATE_DEPTH Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calculate contig depth from BAM files using jgi_summarize_bam_contig_depths
    
    Supports both per-sample and co-assembly modes
    
    Input:
        tuple val(meta), path(contigs), path(bam)
    
    Output:
        depth: tuple val(meta), path(contigs), path(depth_file)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process CALCULATE_DEPTH {
    tag "${meta.id}"
    label 'process_single'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabat2:2.15--h4da6f23_2' :
        'quay.io/biocontainers/metabat2:2.15--h4da6f23_2' }"

    input:
    tuple val(meta), path(contigs), path(bam)

    output:
    tuple val(meta), path(contigs), path("${meta.id}_depth.txt"), emit: depth
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    
    // Input validation
    if (!meta || !meta.id) {
        error "CALCULATE_DEPTH: meta.id is required"
    }
    if (!bam) {
        error "CALCULATE_DEPTH: No BAM file(s) provided for ${meta.id}"
    }
    
    """
    set -euo pipefail
    
    jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam}
    mv depth.txt ${prefix}_depth.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$(jgi_summarize_bam_contig_depths --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """
    
    stub:
    prefix = meta.id
    """
    # Create stub depth file
    cat > ${prefix}_depth.txt << 'EOF'
contigName	contigLen	totalAvgDepth	sample1.bam	sample1.bam-var
contig_1	5000	25.5	25.5	10.2
contig_2	3000	18.3	18.3	8.1
EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: 2.15
    END_VERSIONS
    """
}
