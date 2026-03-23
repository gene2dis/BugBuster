/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METAWRAP Bin Refinement Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Bin refinement using MetaWRAP
    
    Consolidates bins from multiple binners (2 or 3 of: MetaBAT2, SemiBin, COMEBin)
    and refines them based on completeness and contamination thresholds
    
    Input:
        tuple val(meta), path(bin_dirs) - bin_dirs is a list of 2-3 bin directories
    
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

    // Publish refined bins using 'move' mode to free disk space immediately
    publishDir(
        path: "${params.output}/bins/${meta.id}",
        mode: 'move',
        enabled: params.store_refined_bins,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    )

    input:
    tuple val(meta), path(bin_dirs)

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
    
    # Prepare bin directories from each binner and build MetaWRAP flags
    flags=(A B C)
    flag_idx=0
    bin_flags=""
    prep_dirs=""
    
    for bin_dir in ${bin_dirs}; do
        dir_name=\$(basename "\$bin_dir")
        prep_dir="prep_\${flags[\$flag_idx]}"
        mkdir -p "\$prep_dir"
        
        if [[ "\$dir_name" == *metabat* ]]; then
            # Handle MetaBAT2 bins (may be gzipped)
            if [[ -d "\$bin_dir" ]]; then
                for bin in "\$bin_dir"/*.fa* "\$bin_dir"/*.fasta*; do
                    [[ -e "\$bin" ]] || continue
                    if [[ "\$bin" == *.fa.gz ]] || [[ "\$bin" == *.fasta.gz ]]; then
                        gunzip -c "\$bin" > "\$prep_dir/\$(basename "\$bin" .gz)"
                    elif [[ "\$bin" == *.fa ]] || [[ "\$bin" == *.fasta ]]; then
                        cp "\$bin" "\$prep_dir/"
                    fi
                done
            fi
        else
            # SemiBin or COMEBin bins - simple copy
            cp -r "\$bin_dir"/* "\$prep_dir/" 2>/dev/null || true
        fi
        
        bin_flags="\$bin_flags -\${flags[\$flag_idx]} \$prep_dir"
        prep_dirs="\$prep_dirs \$prep_dir"
        flag_idx=\$((flag_idx + 1))
    done
    
    # Count total bins across all prepared directories
    total_bins=\$(find \$prep_dirs -type f \\( -name "*.fa" -o -name "*.fasta" \\) 2>/dev/null | wc -l)
    
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
            \$bin_flags \\
            -c ${completeness} \\
            -x ${contamination}
        
        # Organize output
        chmod 777 Refined_bins
        mv Refined_bins/metawrap_${completeness}_${contamination}_bins ${prefix}_metawrap_${completeness}_${contamination}_bins
    fi
    
    # Cleanup
    rm -rf \$prep_dirs
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metawrap: \$(metawrap --version 2>&1 | head -n1 || echo "1.2")
    END_VERSIONS
    """
    
    stub:
    prefix = meta.id
    completeness = params.metawrap_completeness
    contamination = params.metawrap_contamination
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
