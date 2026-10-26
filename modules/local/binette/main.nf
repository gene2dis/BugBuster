/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BINETTE Bin Refinement Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Bin refinement using Binette

    Builds intermediate bins from the bins of two or three binners (MetaBAT2,
    SemiBin, COMEBin), scores every input and intermediate bin with CheckM2 and
    keeps the best non-overlapping set that reaches params.bin_min_completeness
    and params.bin_max_contamination (Phase 7a benchmark,
    docs/benchmarks/binning_2026.md)

    Input:
        tuple val(meta), path(bin_dirs), path(contigs), path(checkm2_db)

    Output:
        bins: tuple val(meta), path(binette_bins)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process BINETTE {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/binette:1.2.1--pyh7e72e81_0' :
        'quay.io/biocontainers/binette:1.2.1--pyh7e72e81_0' }"

    // Publishing is configured in config/modules.config (withName: 'BINETTE')

    input:
    tuple val(meta), path(bin_dirs), path(contigs), path(checkm2_db)

    output:
    tuple val(meta), path("${meta.id}_binette_bins"), emit: bins
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    def args = task.ext.args ?: ''
    def completeness = task.ext.min_completeness != null ? task.ext.min_completeness : params.bin_min_completeness
    def contamination = task.ext.max_contamination != null ? task.ext.max_contamination : params.bin_max_contamination

    // Input validation
    if (!meta || !meta.id) {
        error "BINETTE: meta.id is required"
    }

    """
    set -euo pipefail

    out_dir=${prefix}_binette_bins

    # Only the binners that found bins go to Binette: a binner directory
    # holding only a note (SKIPPED.txt, FAILED.txt) is left out
    sets=()
    for bin_dir in ${bin_dirs}; do
        if find -L "\$bin_dir" -maxdepth 1 -type f \\( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \\) | grep -q .; then
            sets+=("\$bin_dir")
        fi
    done

    if [[ \${#sets[@]} -eq 0 ]]; then
        echo "WARNING: No bins found from any binner. Skipping Binette refinement for ${prefix}." >&2
        mkdir -p "\$out_dir"
        echo "No bins available for refinement" > "\$out_dir/SKIPPED.txt"
    else
        echo "Found \${#sets[@]} bin sets. Running Binette refinement..."
        binette \\
            --bin_dirs "\${sets[@]}" \\
            --contigs ${contigs} \\
            --outdir binette_out \\
            --threads ${task.cpus} \\
            --checkm2_db ${checkm2_db} \\
            --min_completeness ${completeness} \\
            --max_contamination ${contamination} \\
            --prefix binette \\
            --no-progress \\
            ${args}
        mv binette_out/final_bins "\$out_dir"
        if ! find "\$out_dir" -maxdepth 1 -name '*.fa' | grep -q .; then
            echo "WARNING: No bin reached ${completeness}% completeness and ${contamination}% contamination for ${prefix}." >&2
            echo "No bin reached ${completeness}% completeness and ${contamination}% contamination" > "\$out_dir/SKIPPED.txt"
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$(binette --version 2>&1 | sed 's/^Binette //')
        checkm2: \$(checkm2 --version 2>&1 | sed 's/checkm2: version //')
    END_VERSIONS
    """

    stub:
    prefix = meta.id
    def completeness = task.ext.min_completeness != null ? task.ext.min_completeness : params.bin_min_completeness
    def contamination = task.ext.max_contamination != null ? task.ext.max_contamination : params.bin_max_contamination
    """
    mkdir -p ${prefix}_binette_bins

    # Record the resolved cutoffs so stub-level pipeline tests can pin the
    # parameter wiring (including the deprecated --metawrap_* aliases)
    echo "binette --min_completeness ${completeness} --max_contamination ${contamination}" \\
        > ${prefix}_binette_bins/binette_args.txt

    # Create stub refined bin files
    cat > ${prefix}_binette_bins/binette_bin1.fa << 'EOF'
    >contig_1
    ACGTACGTACGTACGTACGTACGTACGTACGT
    >contig_2
    TGCATGCATGCATGCATGCATGCATGCATGCA
    EOF

    cat > ${prefix}_binette_bins/binette_bin2.fa << 'EOF'
    >contig_3
    GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: 1.2.1
        checkm2: 1.1.0
    END_VERSIONS
    """
}
