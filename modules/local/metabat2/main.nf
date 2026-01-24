/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METABAT2 Binning Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenomic binning using MetaBAT2
    
    Supports both per-sample and co-assembly modes
    
    Input:
        tuple val(meta), path(contigs), path(depth)
    
    Output:
        bins: tuple val(meta), path(metabat_bins)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process METABAT2 {
    tag "${meta.id}"
    label 'process_medium'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabat2:2.15--h4da6f23_2' :
        'quay.io/biocontainers/metabat2:2.15--h4da6f23_2' }"

    input:
    tuple val(meta), path(contigs), path(depth)

    output:
    tuple val(meta), path("${meta.id}_metabat_bins"), emit: bins
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    def args = task.ext.args ?: ''
    
    // Input validation
    if (!meta || !meta.id) {
        error "METABAT2: meta.id is required"
    }
    if (!contigs) {
        error "METABAT2: No contig file provided for sample ${meta.id}"
    }
    if (!depth) {
        error "METABAT2: No depth file provided for sample ${meta.id}"
    }
    
    """
    set -euo pipefail
    
    # Validate input files exist
    if [ ! -f "${contigs}" ]; then
        echo "ERROR: Contig file not found: ${contigs}" >&2
        exit 1
    fi
    
    if [ ! -f "${depth}" ]; then
        echo "ERROR: Depth file not found: ${depth}" >&2
        exit 1
    fi
    
    # Run MetaBAT2
    metabat2 -i ${contigs} \\
             -a ${depth} \\
             -o ${prefix}_metabat_bins/${prefix}_metabat_bin \\
             ${args} \\
             -t ${task.cpus}
    
    # Generate version file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$(metabat2 --help 2>&1 | grep 'version' | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """
    
    stub:
    prefix = meta.id
    """
    mkdir -p ${prefix}_metabat_bins
    
    # Create stub bin files
    cat > ${prefix}_metabat_bins/${prefix}_metabat_bin.1.fa << 'EOF'
>contig_1
ACGTACGTACGTACGTACGTACGTACGTACGT
>contig_2
TGCATGCATGCATGCATGCATGCATGCATGCA
EOF
    
    cat > ${prefix}_metabat_bins/${prefix}_metabat_bin.2.fa << 'EOF'
>contig_3
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: 2.15
    END_VERSIONS
    """
}
