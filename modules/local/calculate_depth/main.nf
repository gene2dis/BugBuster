/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALCULATE_DEPTH Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calculate contig depth from BAM files using jgi_summarize_bam_contig_depths
    (MetaBAT2's input and BIN_DEPTH's MAG depth; runs whenever binning does)
    
    One BAM per sample assembly; under co-assembly one BAM per sample
    (against the co-assembly), giving one depth column per sample
    
    Input:
        tuple val(meta), path(contigs), path(bam) - bam: one BAM or a list
    
    Output:
        depth: tuple val(meta), path(contigs), path(depth_file)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process CALCULATE_DEPTH {
    tag "${meta.id}"
    label 'process_single'
    
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabat2:2.18--h38e344b_2' :
        'quay.io/biocontainers/metabat2:2.18--h38e344b_2' }"

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
        metabat2: \$(jgi_summarize_bam_contig_depths 2>&1 | grep -o 'jgi_summarize_bam_contig_depths [0-9][0-9.]*' | head -n 1 | sed 's/.* //')
    END_VERSIONS
    """
    
    stub:
    prefix = meta.id
    """
    # Stub depth file: jgi's header, one depth and variance column per BAM
    header="contigName\tcontigLen\ttotalAvgDepth"
    for b in ${bam}; do
        header="\${header}\t\$(basename "\$b" .bam)\t\$(basename "\$b" .bam)-var"
    done
    printf "\${header}\n" > ${prefix}_depth.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: 2.18
    END_VERSIONS
    """
}
