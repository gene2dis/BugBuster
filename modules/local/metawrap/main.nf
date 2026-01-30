/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METAWRAP Bin Refinement Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Bin refinement using MetaWRAP
    
    Consolidates bins from multiple binners (MetaBAT2, SemiBin, COMEBin)
    and refines them based on completeness and contamination thresholds
    
    Input:
        tuple val(meta), path(metabat2_bins), path(semibin_bins), path(comebin_bins)
    
    Output:
        bins: tuple val(meta), path(metawrap_bins)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process METAWRAP {
    tag "${meta.id}"
    label 'process_high'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metawrap:1.2--hdfd78af_2' :
        'quay.io/ffuentessantander/metawrap:1.2' }"

    // Use storeDir for refined bins to enable immediate work dir cleanup
    // This stores final bins permanently and cleans work directory automatically
    storeDir params.store_refined_bins ? "${params.output}/bins/${meta.id}" : null

    input:
    tuple val(meta), path(metabat2_bins), path(semibin_bins), path(comebin_bins)

    output:
    tuple val(meta), path("${meta.id}_metawrap_${params.metawrap_completeness}_${params.metawrap_contamination}_bins"), emit: bins
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    completeness = params.metawrap_completeness
    contamination = params.metawrap_contamination
    
    // Input validation
    if (!meta || !meta.id) {
        error "METAWRAP: meta.id is required"
    }
    
    """
    set -euo pipefail
    
    # Prepare bin directories
    mkdir -p metabat_wp_bins semibin_wp_bins comebin_wp_bins
    
    # Handle MetaBAT2 bins (may be gzipped)
    if [[ -d "${metabat2_bins}" ]]; then
        for bin in ${metabat2_bins}/*.fa* ${metabat2_bins}/*.fasta*; do
            [[ -e "\$bin" ]] || continue
            if [[ "\$bin" == *.fa.gz ]] || [[ "\$bin" == *.fasta.gz ]]; then
                gunzip -c "\$bin" > "metabat_wp_bins/\$(basename "\$bin" .gz)"
            elif [[ "\$bin" == *.fa ]] || [[ "\$bin" == *.fasta ]]; then
                cp "\$bin" metabat_wp_bins/
            fi
        done
    else
        # Handle single file input
        if [[ "${metabat2_bins}" == *.fa.gz ]] || [[ "${metabat2_bins}" == *.fasta.gz ]]; then
            gunzip -c "${metabat2_bins}" > "metabat_wp_bins/\$(basename "${metabat2_bins}" .gz)"
        elif [[ "${metabat2_bins}" == *.fa ]] || [[ "${metabat2_bins}" == *.fasta ]]; then
            cp "${metabat2_bins}" metabat_wp_bins/
        fi
    fi
    
    # Copy SemiBin and COMEBin bins (avoid dereferencing to prevent infinite loops)
    cp -r ${semibin_bins}/* semibin_wp_bins/ 2>/dev/null || true
    cp -r ${comebin_bins}/* comebin_wp_bins/ 2>/dev/null || true
    
    # Count total bins across all binners
    total_bins=\$(find metabat_wp_bins semibin_wp_bins comebin_wp_bins -type f \\( -name "*.fa" -o -name "*.fasta" \\) 2>/dev/null | wc -l)
    
    if [[ \$total_bins -eq 0 ]]; then
        echo "WARNING: No bins found from any binner. Skipping MetaWRAP refinement."
        echo "Creating empty output directory for ${prefix}"
        mkdir -p ${prefix}_metawrap_${completeness}_${contamination}_bins
        echo "No bins available for refinement" > ${prefix}_metawrap_${completeness}_${contamination}_bins/SKIPPED.txt
    else
        echo "Found \$total_bins bins total. Running MetaWRAP refinement..."
        
        # Run MetaWRAP bin refinement
        metawrap bin_refinement \\
            -o Refined_bins \\
            -t ${task.cpus} \\
            -A metabat_wp_bins \\
            -B semibin_wp_bins \\
            -C comebin_wp_bins \\
            -c ${completeness} \\
            -x ${contamination}
        
        # Organize output
        chmod 777 Refined_bins
        mv Refined_bins/metawrap_${completeness}_${contamination}_bins ${prefix}_metawrap_${completeness}_${contamination}_bins
    fi
    
    # Cleanup
    rm -r metabat_wp_bins semibin_wp_bins comebin_wp_bins
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metawrap: \$(metawrap --version 2>&1 | head -n1 || echo "1.2")
    END_VERSIONS
    """
    
    stub:
    prefix = meta.id
    """
    mkdir -p ${prefix}_metawrap_${completeness}_${contamination}_bins
    
    # Create stub refined bin files
    cat > ${prefix}_metawrap_${completeness}_${contamination}_bins/bin.1.fa << 'EOF'
>contig_1
ACGTACGTACGTACGTACGTACGTACGTACGT
>contig_2
TGCATGCATGCATGCATGCATGCATGCATGCA
EOF
    
    cat > ${prefix}_metawrap_${completeness}_${contamination}_bins/bin.2.fa << 'EOF'
>contig_3
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metawrap: 1.2
    END_VERSIONS
    """
}
