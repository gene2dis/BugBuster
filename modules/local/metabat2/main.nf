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
    
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabat2:2.18--h38e344b_2' :
        'quay.io/biocontainers/metabat2:2.18--h38e344b_2' }"

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
    
    # Run MetaBAT2. With contigs but none of --metabat_minContig bp it stops
    # with exit 1 ("There were no large target contigs"), which ended the run
    # for a sample that just has no bins: write the bin directory with a
    # SKIPPED.txt instead, as SEMIBIN and COMEBIN do (Binette, CHECKM2_BATCH
    # and GTDB_TK_BATCH skip a directory without .fa files; roadmap Phase 21,
    # RV10). Every other failure still fails the task
    status=0
    metabat2 -i ${contigs} \\
             -a ${depth} \\
             -o ${prefix}_metabat_bins/${prefix}_metabat_bin \\
             ${args} \\
             -t ${task.cpus} > metabat2.log 2>&1 || status=\$?
    cat metabat2.log >&2
    if [ "\$status" -ne 0 ]; then
        if grep -q 'There were no large target contigs' metabat2.log; then
            echo "WARNING: no contig of ${prefix} reaches MetaBAT2's minimum contig length; MetaBAT2 skipped" >&2
            mkdir -p ${prefix}_metabat_bins
            echo "No contig reaches MetaBAT2's minimum contig length (--metabat_minContig) - MetaBAT2 skipped" > ${prefix}_metabat_bins/SKIPPED.txt
        else
            exit "\$status"
        fi
    fi
    rm -f metabat2.log

    # MetaBAT2 2.18 writes depth fields after the contig name in bin FASTA
    # headers; refiners and the bin depth step match contigs by the whole
    # header, so keep the name only (2.15 wrote the bare name)
    find ${prefix}_metabat_bins -name '*.fa' -exec sed -i '/^>/s/[[:space:]].*//' {} +
    
    # Generate version file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$(metabat2 --help 2>&1 | grep -o 'version [0-9:]*[0-9]\\.[0-9.]*' | head -n 1 | sed 's/.*[ :]//')
    END_VERSIONS
    """
    
    stub:
    prefix = meta.id
    def args = task.ext.args ?: ''
    """
    mkdir -p ${prefix}_metabat_bins

    # Record the resolved ext.args so stub-level pipeline tests can pin the
    # config wiring (audit #14: conflicting selectors dropped params silently)
    echo "metabat2 ${args}" > ${prefix}_metabat_bins/metabat2_args.txt

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
        metabat2: 2.18
    END_VERSIONS
    """
}
