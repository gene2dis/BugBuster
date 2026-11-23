/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GTDB_TK_BATCH Taxonomic Classification Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Batched taxonomic classification of bins using GTDB-Tk
    
    Processes all bins from all samples in a single GTDB-Tk run for improved performance.
    Database is loaded once instead of N times (N = number of samples).
    Removes maxForks limitation for better parallelization.
    
    Input:
        tuple val(meta_list), path(all_bins, stageAs: 'bin_?/*'), path(gtdbtk_db)
    
    Output:
        gtdb_tk: tuple val(meta_list), path(gtdbtk_files) - per sample
        report: tuple val(meta_list), path(gtdbtk_files) - per sample
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process GTDB_TK_BATCH {
    tag "batch_all_samples"
    label 'process_high'
    label 'process_high_memory'
    
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.8.0--pyhdfd78af_0' :
        'quay.io/biocontainers/gtdbtk:2.8.0--pyhdfd78af_0' }"

    input:
    tuple val(meta_list), path(all_bins, stageAs: 'bin_?/*'), path(gtdbtk_db)

    output:
    tuple val(meta_list), path("*_gtdbtk_*"), emit: gtdb_tk
    tuple val(meta_list), path("*_gtdbtk_*"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def meta_ids = meta_list.collect { m -> m instanceof Map ? (m.id ?: m['id']) : m.toString() }.unique().join(' ')
    // pplacer threads (performance review M10): 4 placed genomes 1.5-2x
    // faster than 1 for ~7 % more memory; never more than the task's CPUs
    def pplacer_cpus = Math.min(params.gtdbtk_pplacer_cpus as Integer, task.cpus as Integer)
    // The sample of each staged bin directory, from the meta list (one meta
    // per path, in the order Nextflow stages them as bin_1, bin_2, ...), so no
    // id is recovered from a directory name or the resolved work path
    // (roadmap Phase 20, RV5)
    def bin_paths = all_bins instanceof List ? all_bins : [all_bins]
    if (bin_paths.size() != meta_list.size()) {
        error "GTDB_TK_BATCH: ${bin_paths.size()} bin directories for ${meta_list.size()} metas"
    }
    // one line of quoted pairs: a multi-line value would break the script's
    // indentation stripping (R4-04)
    def staged_ids = [bin_paths, meta_list].transpose().collect { p, m -> "'${p}' '${m instanceof Map ? m.id : m}'" }.join(' ')
    """
    set -euo pipefail

    # staged bin directory -> sample id
    printf '%s\\t%s\\n' ${staged_ids} > staged_ids.tsv

    # Create combined directory with prefixed bin names
    mkdir -p combined_bins

    total_bins=0

    # Process ALL staged bin directories; the sample id comes from staged_ids.tsv
    # Bins are staged as bin_1/*, bin_2/*, etc. by Nextflow. The input channel
    # already carries exactly the dirs to classify (single-binner outputs, or
    # Binette-refined dirs), so every dir must be recognized; identify_bin_dir.sh
    # fails loudly on naming drift (audit #6).
    for staged_dir in bin_*/; do
        staged_dir=\${staged_dir%/}  # Remove trailing slash

        # Get the actual bin directory inside
        bin_dir=\$(find "\$staged_dir" -mindepth 1 -maxdepth 1 -type d -o -type l | head -1)
        [ -z "\$bin_dir" ] && continue

        staged_id=\$(awk -F '\\t' -v d="\$bin_dir" '\$1 == d {print \$2}' staged_ids.tsv)
        [ -n "\$staged_id" ] || { echo "ERROR: \$bin_dir is not in staged_ids.tsv (staged inputs and the meta list disagree)" >&2; exit 1; }
        id_line=\$(identify_bin_dir.sh "\$bin_dir" "\$staged_id")
        sample_id=\${id_line%%\$'\\t'*}

        echo "Processing \$(basename "\$bin_dir") -> sample_id=\$sample_id"

        # Find all .fa/.fasta/.fna files recursively
        bin_count=\$(find -L "\$bin_dir" -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) 2>/dev/null | wc -l)

        if [ "\$bin_count" -gt 0 ]; then
            echo "Found \$bin_count bins for sample \${sample_id}"
            # Copy bins with prefixed names to track origin
            find -L "\$bin_dir" -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) | while read bin_file; do
                bin_name=\$(basename "\$bin_file")
                # named .fa whatever the extension: GTDB-Tk runs with
                # --extension .fa, so a .fasta / .fna bin was counted but not
                # classified (R4-12)
                cp "\$bin_file" "combined_bins/\${sample_id}__\${bin_name%.*}.fa"
            done
            total_bins=\$((total_bins + bin_count))
        fi
    done
    
    if [ "\$total_bins" -eq 0 ]; then
        echo "WARNING: No genome bins found across all samples. Creating empty reports."
        
        # Create empty reports for each sample
        for meta_id in ${meta_ids}; do
            echo -e "user_genome\\tclassification\\tclosest_genome_reference\\tclosest_genome_reference_radius\\tclosest_genome_taxonomy\\tclosest_genome_ani\\tclosest_genome_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references(genome_id,species_name,radius,ANI,AF)\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings" > \${meta_id}_gtdbtk_bac120.tsv
            echo -e "user_genome\\tclassification\\tclosest_genome_reference\\tclosest_genome_reference_radius\\tclosest_genome_taxonomy\\tclosest_genome_ani\\tclosest_genome_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references(genome_id,species_name,radius,ANI,AF)\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings" > \${meta_id}_gtdbtk_ar53.tsv
        done
    else
        echo "Found \$total_bins genome bins across all samples. Running GTDB-Tk classification..."
        
        export GTDBTK_DATA_PATH="${gtdbtk_db}"

        # GTDB-Tk 2.7 removed --skip_ani_screen/--mash_db: the ANI screen
        # always runs, against the skani DB bundled in the R232 data package
        gtdbtk classify_wf \\
            --genome_dir combined_bins \\
            --out_dir gtdbtk_output \\
            --cpus ${task.cpus} \\
            --extension .fa \\
            --pplacer_cpus ${pplacer_cpus}
        
        # Split the combined reports by sample prefix (<sample>__<bin>):
        # bin/gtdbtk_split.py stops on a bin whose prefix is not a sample of
        # the batch instead of dropping its row (roadmap Phase 20, RV3), and
        # writes header-only tables for a domain GTDB-Tk wrote no summary for
        # (no genome placed in it): ar53 for bacteria-only bins, bac120 for
        # archaea-only bins (roadmap Phase 21, RV9; review R4-03)
        if [ -f gtdbtk_output/classify/gtdbtk.bac120.summary.tsv ] || [ -f gtdbtk_output/classify/gtdbtk.ar53.summary.tsv ]; then
            gtdbtk_split.py --samples ${meta_ids} --domain bac120 --summary gtdbtk_output/classify/gtdbtk.bac120.summary.tsv
            gtdbtk_split.py --samples ${meta_ids} --domain ar53 --summary gtdbtk_output/classify/gtdbtk.ar53.summary.tsv
        else
            echo "ERROR: GTDB-Tk wrote neither gtdbtk.bac120.summary.tsv nor gtdbtk.ar53.summary.tsv for \$total_bins bins; the end of its log:" >&2
            tail -n 50 gtdbtk_output/gtdbtk.log >&2 || true
            exit 1
        fi
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(gtdbtk --version 2>&1 | head -n 1 | sed 's/gtdbtk: version //; s/ Copyright.*//')
    END_VERSIONS
    """
    
    stub:
    def meta_ids = meta_list.collect { m -> m instanceof Map ? (m.id ?: m['id']) : m.toString() }.unique().join(' ')
    """
    # Create stub taxonomy files for each sample
    for meta_id in ${meta_ids}; do
        cat > \${meta_id}_gtdbtk_bac120.tsv << 'EOF'
    user_genome	classification
    bin.1	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
    bin.2	d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis
    EOF
        
        touch \${meta_id}_gtdbtk_ar53.tsv
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: 2.8.0
    END_VERSIONS
    """
}
