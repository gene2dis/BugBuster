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
    
    # Check if bins directory contains actual genome files
    bin_count=\$(find -L ${metawrap} -maxdepth 1 -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) 2>/dev/null | wc -l)
    
    if [ "\$bin_count" -eq 0 ]; then
        echo "WARNING: No genome bins found in ${metawrap}. Skipping GTDB-Tk classification."
        
        # Create empty output files with headers
        echo -e "user_genome\\tclassification\\tfastani_reference\\tfastani_reference_radius\\tfastani_taxonomy\\tfastani_ani\\tfastani_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references(genome_id,species_name,radius,ANI,AF)\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings" > ${prefix}_gtdbtk_bac120.tsv
        echo -e "user_genome\\tclassification\\tfastani_reference\\tfastani_reference_radius\\tfastani_taxonomy\\tfastani_ani\\tfastani_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references(genome_id,species_name,radius,ANI,AF)\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings" > ${prefix}_gtdbtk_ar53.tsv
    else
        echo "Found \$bin_count genome bins. Running GTDB-Tk classification..."
        
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
