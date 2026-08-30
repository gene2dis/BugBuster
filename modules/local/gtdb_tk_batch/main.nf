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
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.7.2--pyhdfd78af_1' :
        'quay.io/biocontainers/gtdbtk:2.7.2--pyhdfd78af_1' }"

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
    def metawrap_suffix = task.ext.metawrap_dir_suffix ?: "metawrap_${params.metawrap_completeness}_${params.metawrap_contamination}_bins"
    """
    set -euo pipefail

    # Create combined directory with prefixed bin names
    mkdir -p combined_bins

    total_bins=0

    # Process ALL staged bin directories - parse sample_id from directory names
    # Bins are staged as bin_1/*, bin_2/*, etc. by Nextflow. The input channel
    # already carries exactly the dirs to classify (single-binner outputs, or
    # MetaWRAP-refined dirs), so every dir must be recognized; identify_bin_dir.sh
    # fails loudly on naming drift (audit #6).
    for staged_dir in bin_*/; do
        staged_dir=\${staged_dir%/}  # Remove trailing slash

        # Get the actual bin directory inside
        bin_dir=\$(find "\$staged_dir" -mindepth 1 -maxdepth 1 -type d -o -type l | head -1)
        [ -z "\$bin_dir" ] && continue

        id_line=\$(identify_bin_dir.sh "\$bin_dir" "${metawrap_suffix}")
        sample_id=\${id_line%%\$'\\t'*}

        echo "Processing \$(basename "\$bin_dir") -> sample_id=\$sample_id"

        # Find all .fa/.fasta/.fna files recursively
        bin_count=\$(find -L "\$bin_dir" -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) 2>/dev/null | wc -l)

        if [ "\$bin_count" -gt 0 ]; then
            echo "Found \$bin_count bins for sample \${sample_id}"
            # Copy bins with prefixed names to track origin
            find -L "\$bin_dir" -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) | while read bin_file; do
                bin_name=\$(basename "\$bin_file")
                cp "\$bin_file" "combined_bins/\${sample_id}__\${bin_name}"
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
            --pplacer_cpus 1
        
        # Split the combined reports by sample prefix
        if [ -f gtdbtk_output/classify/gtdbtk.bac120.summary.tsv ]; then
            python3 <<-PYEOF
			import sys
			import os
			
			meta_ids = "${meta_ids}".split()
			
			# Process bacteria summary
			if os.path.exists("gtdbtk_output/classify/gtdbtk.bac120.summary.tsv"):
			    with open("gtdbtk_output/classify/gtdbtk.bac120.summary.tsv", 'r') as f:
			        header = f.readline()
			        lines = f.readlines()
			    
			    # Create per-sample reports
			    sample_reports = {meta_id: [] for meta_id in meta_ids}
			    
			    for line in lines:
			        if line.strip():
			            # Extract sample ID from bin name (format: sampleID__binname.fa)
			            bin_name = line.split('\t')[0]
			            if '__' in bin_name:
			                sample_id = bin_name.split('__')[0]
			                if sample_id in sample_reports:
			                    # Remove the sample prefix from bin name in output
			                    original_bin_name = '__'.join(bin_name.split('__')[1:])
			                    modified_line = line.replace(bin_name, original_bin_name, 1)
			                    sample_reports[sample_id].append(modified_line)
			    
			    # Write per-sample reports
			    for meta_id in meta_ids:
			        with open(f"{meta_id}_gtdbtk_bac120.tsv", 'w') as f:
			            f.write(header)
			            if sample_reports[meta_id]:
			                f.writelines(sample_reports[meta_id])
			
			# Process archaea summary
			if os.path.exists("gtdbtk_output/classify/gtdbtk.ar53.summary.tsv"):
			    with open("gtdbtk_output/classify/gtdbtk.ar53.summary.tsv", 'r') as f:
			        header = f.readline()
			        lines = f.readlines()
			    
			    # Create per-sample reports
			    sample_reports = {meta_id: [] for meta_id in meta_ids}
			    
			    for line in lines:
			        if line.strip():
			            # Extract sample ID from bin name (format: sampleID__binname.fa)
			            bin_name = line.split('\t')[0]
			            if '__' in bin_name:
			                sample_id = bin_name.split('__')[0]
			                if sample_id in sample_reports:
			                    # Remove the sample prefix from bin name in output
			                    original_bin_name = '__'.join(bin_name.split('__')[1:])
			                    modified_line = line.replace(bin_name, original_bin_name, 1)
			                    sample_reports[sample_id].append(modified_line)
			    
			    # Write per-sample reports
			    for meta_id in meta_ids:
			        with open(f"{meta_id}_gtdbtk_ar53.tsv", 'w') as f:
			            f.write(header)
			            if sample_reports[meta_id]:
			                f.writelines(sample_reports[meta_id])
			else:
			    # Create empty ar53 files if no archaea found
			    for meta_id in meta_ids:
			        if not os.path.exists(f"{meta_id}_gtdbtk_ar53.tsv"):
			            with open(f"{meta_id}_gtdbtk_ar53.tsv", 'w') as f:
			                f.write("user_genome\\tclassification\\tclosest_genome_reference\\tclosest_genome_reference_radius\\tclosest_genome_taxonomy\\tclosest_genome_ani\\tclosest_genome_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references(genome_id,species_name,radius,ANI,AF)\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings\\n")
			PYEOF
        fi
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(gtdbtk --version 2>&1 | head -n 1 | sed 's/gtdbtk: version //; s/ Copyright.*//')
    END_VERSIONS
    """
    
    stub:
    def meta_ids = meta_list.collect { m -> m instanceof Map ? (m.id ?: m['id']) : m.toString() }.join(' ')
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
        gtdbtk: 2.7.2
    END_VERSIONS
    """
}
