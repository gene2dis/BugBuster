/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GTDB-TK Taxonomic Classification Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Taxonomic classification of bins using GTDB-Tk
    
    Supports both per-sample and co-assembly modes
    
    Input:
        tuple val(meta), path(metawrap), path(gtdbtk_db)
    
    Output:
        gtdb_tk: tuple val(meta), path(gtdbtk_files)
        report: tuple val(meta), path(gtdbtk_files)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process GTDB_TK {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'
    maxForks 2
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.5.2--pyh1f0d9b5_0' :
        'quay.io/biocontainers/gtdbtk:2.5.2--pyh1f0d9b5_0' }"

    input:
    tuple val(meta), path(metawrap), path(gtdbtk_db)

    output:
    tuple val(meta), path("*_gtdbtk_*"), emit: gtdb_tk
    tuple val(meta), path("*_gtdbtk_*"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    
    // Input validation
    if (!meta || !meta.id) {
        error "GTDB_TK: meta.id is required"
    }
    if (!metawrap) {
        error "GTDB_TK: No bin directory provided for ${meta.id}"
    }
    if (!gtdbtk_db) {
        error "GTDB_TK: GTDB-Tk database path is required"
    }
    
    """
    set -euo pipefail
    
    export GTDBTK_DATA_PATH="${gtdbtk_db}"
    
    gtdbtk classify_wf \\
        --genome_dir ${metawrap} \\
        --out_dir ${prefix}_gtdb \\
        --cpus ${task.cpus} \\
        --skip_ani_screen \\
        --extension .fa \\
        --pplacer_cpus 1
    
    mv ${prefix}_gtdb/classify/gtdbtk.bac120.summary.tsv ${prefix}_gtdbtk_bac120.tsv
    
    if [ -f ${prefix}_gtdb/classify/gtdbtk.ar53.summary.tsv ]; then
        mv ${prefix}_gtdb/classify/gtdbtk.ar53.summary.tsv ${prefix}_gtdbtk_ar53.tsv
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(gtdbtk --version 2>&1 | sed 's/gtdbtk: version //')
    END_VERSIONS
    """
    
    stub:
    prefix = meta.id
    """
    # Create stub taxonomy files
    cat > ${prefix}_gtdbtk_bac120.tsv << 'EOF'
user_genome	classification
bin.1	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
bin.2	d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis
EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: 2.5.2
    END_VERSIONS
    """
}
