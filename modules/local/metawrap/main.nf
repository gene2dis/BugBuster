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

    // Publishing is configured in config/modules.config (withName: 'METAWRAP')

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

    out_dir=${prefix}_metawrap_${completeness}_${contamination}_bins
    final=metawrap_${completeness}_${contamination}_bins

    # One prep dir per binner holding only its FASTA bins (MetaBAT2 bins may
    # be gzipped). A note such as SKIPPED.txt in a binner's empty directory is
    # left out: MetaWRAP would count it as a bin
    n=0
    for bin_dir in ${bin_dirs}; do
        n=\$((n + 1))
        mkdir -p prep_\$n
        for bin in "\$bin_dir"/*.fa "\$bin_dir"/*.fasta "\$bin_dir"/*.fa.gz "\$bin_dir"/*.fasta.gz; do
            [[ -e "\$bin" ]] || continue
            case "\$bin" in
                *.gz) gunzip -c "\$bin" > "prep_\$n/\$(basename "\$bin" .gz)" ;;
                *)    cp "\$bin" "prep_\$n/" ;;
            esac
        done
    done

    # MetaWRAP stops on an empty bin set: pass only the binners that found
    # bins, lettered A, B, C in input order
    letters=(A B C)
    bin_flags=""
    n_sets=0
    for prep_dir in prep_*; do
        [[ -n "\$(ls -A "\$prep_dir")" ]] || continue
        bin_flags="\$bin_flags -\${letters[\$n_sets]} \$prep_dir"
        n_sets=\$((n_sets + 1))
    done

    # MetaWRAP exits 0 even when it stops on an error, so its outcome is read
    # from the directories it leaves
    refine() {
        metawrap bin_refinement \\
            -o "\$1" \\
            -t ${task.cpus} \\
            \$bin_flags \\
            -c ${completeness} \\
            -x ${contamination} \\
            \${2:-}
        chmod 777 "\$1"
    }
    has_empty_dir() {
        local dir
        for dir in "\$@"; do
            if [[ -d "\$dir" && -z "\$(ls -A "\$dir")" ]]; then return 0; fi
        done
        return 1
    }

    if [[ \$n_sets -eq 0 ]]; then
        echo "WARNING: No bins found from any binner. Skipping MetaWRAP refinement for ${prefix}." >&2
        mkdir -p "\$out_dir"
        echo "No bins available for refinement" > "\$out_dir/SKIPPED.txt"
    elif [[ \$n_sets -eq 1 ]]; then
        # MetaWRAP 1.2 cannot refine a single bin set: it stops after
        # CheckM, looking for its consolidated set. Keep the bins that pass
        # its own quality test on its CheckM table (binsA.stats), so refined
        # bins always meet the thresholds
        echo "WARNING: Only one binner found bins for ${prefix}; keeping its bins that reach ${completeness}% completeness and ${contamination}% contamination." >&2
        refine Refined_bins
        if [[ ! -s Refined_bins/binsA.stats ]]; then
            echo "ERROR: MetaWRAP's CheckM run failed for ${prefix}; see its log above." >&2
            exit 1
        fi
        mkdir -p "\$out_dir"
        awk -v c=${completeness} -v x=${contamination} 'NR > 1 && \$2 >= c && \$2 <= 100 && \$3 >= 0 && \$3 <= x { print \$1 }' Refined_bins/binsA.stats \\
            | while read -r bin; do cp "Refined_bins/binsA/\$bin.fa" "\$out_dir/"; done
        if [[ -z "\$(ls -A "\$out_dir")" ]]; then
            echo "No bin reached ${completeness}% completeness and ${contamination}% contamination" > "\$out_dir/SKIPPED.txt"
        fi
    else
        echo "Found \$n_sets bin sets. Running MetaWRAP refinement..."
        result=Refined_bins
        refine \$result
        # A hybrid set is empty when two binners' bins share too little
        # sequence; MetaWRAP stops on it. Refine the binners' own sets instead
        if [[ ! -d \$result/\$final ]] && has_empty_dir \$result/binsAB \$result/binsAC \$result/binsBC \$result/binsABC; then
            echo "WARNING: An empty hybrid bin set stopped MetaWRAP for ${prefix}; refining the binners' own bin sets only (--skip-refinement)." >&2
            result=Refined_bins_original_sets
            refine \$result --skip-refinement
        fi

        if [[ -d \$result/\$final ]]; then
            mv \$result/\$final "\$out_dir"
        elif has_empty_dir \$result/binsO; then
            # MetaWRAP stops when no bin reaches the thresholds
            echo "WARNING: No bin reached ${completeness}% completeness and ${contamination}% contamination for ${prefix}." >&2
            mkdir -p "\$out_dir"
            echo "No bin reached ${completeness}% completeness and ${contamination}% contamination" > "\$out_dir/SKIPPED.txt"
        else
            echo "ERROR: MetaWRAP bin_refinement failed for ${prefix}; see its log above." >&2
            exit 1
        fi
    fi

    # Cleanup
    rm -rf prep_*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metawrap: \$(metawrap --version 2>&1 | grep -o 'v=[0-9.]*' | sed 's/v=//')
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
