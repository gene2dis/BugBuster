/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BINETTE Bin Refinement Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Bin refinement using Binette

    Builds intermediate bins from the bins of two or three binners (MetaBAT2,
    SemiBin, COMEBin), scores every input and intermediate bin with CheckM2 and
    keeps the best non-overlapping set that reaches params.bin_min_completeness
    and params.bin_max_contamination (Phase 7a benchmark,
    docs/benchmarks/binning_2026.md). Its CheckM2 scores of the input and
    refined bins are written in CHECKM2_BATCH's per-binner report format
    (bin/binette_quality_reports.py), so no separate CheckM2 run is needed
    when bins are refined (performance review M7)

    Input:
        tuple val(meta), path(bin_dirs), path(contigs), path(checkm2_db)

    Output:
        bins: tuple val(meta), path(binette_bins)
        quality_reports: tuple val(meta), path(<id>_<binner>_quality_report.tsv), one per
            staged binner directory plus <id>_binette_quality_report.tsv
        refined_report: tuple val(meta), path(<id>_binette_quality_report.tsv)
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
    tuple val(meta), path("*_quality_report.tsv"), emit: quality_reports
    tuple val(meta), path("${meta.id}_binette_quality_report.tsv"), emit: refined_report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    def args = task.ext.args ?: ''
    // --bin_min_completeness / --bin_max_contamination, else their deprecated
    // --metawrap_* aliases (main.nf stops the run when both are given with
    // different values), else 50 / 10
    def completeness = [params.bin_min_completeness, params.metawrap_completeness, 50].find { v -> v != null }
    def contamination = [params.bin_max_contamination, params.metawrap_contamination, 10].find { v -> v != null }

    // Input validation
    if (!meta || !meta.id) {
        error "BINETTE: meta.id is required"
    }

    """
    set -euo pipefail

    out_dir=${prefix}_binette_bins

    # Directory holds at least one FASTA bin. find stops at the first match:
    # 'find | grep -q' let find die of SIGPIPE under pipefail on a directory
    # with many bins, which read as "no bins"
    has_bins() {
        [[ -n "\$(find -L "\$1" -maxdepth 1 -type f \\( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \\) -print -quit)" ]]
    }

    # Each binner's directory is staged as sets/<binner> (metabat, semibin,
    # comebin): Binette names bin sets after their paths, and sample ids with
    # dots made two binners' sets collide. Only the binners that found bins go
    # to Binette: a directory holding only a note (SKIPPED.txt, FAILED.txt) is
    # left out
    mkdir -p sets
    for bin_dir in ${bin_dirs}; do
        id_line=\$(identify_bin_dir.sh "\$bin_dir")
        ln -s "../\$bin_dir" "sets/\${id_line##*\$'\\t'}"
    done
    with_bins=()
    for set_dir in sets/*; do
        if has_bins "\$set_dir"; then
            with_bins+=("\$set_dir")
        fi
    done

    if [[ \${#with_bins[@]} -eq 0 ]]; then
        echo "WARNING: No bins found from any binner. Skipping Binette refinement for ${prefix}." >&2
        mkdir -p "\$out_dir"
        echo "No bins available for refinement" > "\$out_dir/SKIPPED.txt"
    else
        echo "Found \${#with_bins[@]} bin sets. Running Binette refinement..."
        binette \\
            --bin_dirs "\${with_bins[@]}" \\
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
        if ! has_bins "\$out_dir"; then
            echo "WARNING: No bin reached ${completeness}% completeness and ${contamination}% contamination for ${prefix}." >&2
            echo "No bin reached ${completeness}% completeness and ${contamination}% contamination" > "\$out_dir/SKIPPED.txt"
        fi
    fi

    # CheckM2 reports of every binner's bins and of the refined bins, from
    # Binette's own CheckM2 run (header-only when Binette did not run)
    binette_out_arg=""
    [[ -d binette_out ]] && binette_out_arg="--binette-out binette_out"
    binette_quality_reports.py \\
        --sample ${prefix} \\
        \$binette_out_arg \\
        --final-bins "\$out_dir" \\
        sets/*
    rm -rf binette_out/temporary_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$(binette --version 2>&1 | sed 's/^Binette //')
        checkm2: \$(checkm2 --version 2>&1 | sed 's/checkm2: version //')
    END_VERSIONS
    """

    stub:
    prefix = meta.id
    // --bin_min_completeness / --bin_max_contamination, else their deprecated
    // --metawrap_* aliases (main.nf stops the run when both are given with
    // different values), else 50 / 10
    def completeness = [params.bin_min_completeness, params.metawrap_completeness, 50].find { v -> v != null }
    def contamination = [params.bin_max_contamination, params.metawrap_contamination, 10].find { v -> v != null }
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

    # Header-only CheckM2 reports, one per binner
    mkdir -p sets
    for bin_dir in ${bin_dirs}; do
        id_line=\$(identify_bin_dir.sh "\$bin_dir")
        ln -s "../\$bin_dir" "sets/\${id_line##*\$'\\t'}"
    done
    binette_quality_reports.py --sample ${prefix} --final-bins ${prefix}_binette_bins sets/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: 1.2.1
        checkm2: 1.1.0
    END_VERSIONS
    """
}
