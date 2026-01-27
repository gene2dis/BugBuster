/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECKM2_BATCH Quality Assessment Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Batched bin quality assessment using CheckM2
    
    Processes all bins from all samples in a single CheckM2 run for improved performance.
    Database is loaded once instead of N times (N = number of samples).
    
    Input:
        tuple val(meta_list), path(all_bins, stageAs: 'input_bins/*'), path(checkm_db)
    
    Output:
        all_reports: tuple val(meta), path(quality_reports) - per sample
        metawrap_report: tuple val(meta), path(metawrap_quality_report) - per sample
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process CHECKM2_BATCH {
    tag "batch_all_samples"
    label 'process_medium'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.1.0--pyh7e72e81_1' :
        'quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1' }"

    input:
    tuple val(meta_list), path(all_bins, stageAs: 'bin_?/*'), path(checkm_db)

    output:
    tuple val(meta_list), path("*_*_quality_report.tsv"), emit: all_reports
    tuple val(meta_list), path("*_metawrap_quality_report.tsv"), emit: metawrap_report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def meta_ids = meta_list.collect { m -> m instanceof Map ? (m.id ?: m['id']) : m.toString() }.unique().join(' ')
    """
    set -euo pipefail
    
    # Create organized directory structure
    mkdir -p batch_bins
    
    # Process ALL staged bin directories - parse sample_id and binner from directory names
    # Bins are staged as bin_1/*, bin_2/*, etc. by Nextflow
    for staged_dir in bin_*/; do
        staged_dir=\${staged_dir%/}  # Remove trailing slash
        
        # Get the actual bin directory inside
        bin_dir=\$(find "\$staged_dir" -mindepth 1 -maxdepth 1 -type d -o -type l | head -1)
        [ -z "\$bin_dir" ] && continue
        
        bin_dir_name=\$(basename "\$bin_dir")
        
        # Parse directory name to extract sample_id and binner
        # Pattern: {sample_id}_{binner}_bins or {sample_id}_{binner}_output_bins
        if [[ "\$bin_dir_name" =~ ^(.+)_metabat_bins\$ ]]; then
            sample_id="\${BASH_REMATCH[1]}"
            binner="metabat"
        elif [[ "\$bin_dir_name" =~ ^(.+)_semibin_output_bins\$ ]]; then
            sample_id="\${BASH_REMATCH[1]}"
            binner="semibin"
        elif [[ "\$bin_dir_name" =~ ^(.+)_comebin_bins\$ ]] || [[ "\$bin_dir_name" =~ comebin_res_bins\$ ]]; then
            # COMEBin has nested structure - extract sample from parent or use directory listing
            if [[ "\$bin_dir_name" =~ ^(.+)_comebin ]]; then
                sample_id="\${BASH_REMATCH[1]}"
            else
                # For comebin_res_bins symlink, get sample from symlink target
                sample_id=\$(readlink -f "\$bin_dir" | grep -oP '[^/]+(?=_comebin)' || echo "unknown")
            fi
            binner="comebin"
        elif [[ "\$bin_dir_name" =~ ^(.+)_metawrap_50_10_bins\$ ]]; then
            sample_id="\${BASH_REMATCH[1]}"
            binner="metawrap"
        else
            echo "Skipping unrecognized directory: \$bin_dir_name"
            continue
        fi
        
        echo "Processing \$bin_dir -> sample_id=\$sample_id, binner=\$binner"
        
        mkdir -p "batch_bins/\${sample_id}_\${binner}"
        
        # Find and copy bins with prefixed names
        bin_count=\$(find -L "\$bin_dir" -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) 2>/dev/null | wc -l)
        
        if [ "\$bin_count" -gt 0 ]; then
            echo "Found \$bin_count bins for \${sample_id} \${binner}"
            find -L "\$bin_dir" -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) | while read bin_file; do
                bin_name=\$(basename "\$bin_file")
                cp "\$bin_file" "batch_bins/\${sample_id}_\${binner}/\${sample_id}__\${bin_name}"
            done
        fi
    done
    
    # Run CheckM2 on each binner's combined bins across all samples
    for binner in metabat semibin comebin metawrap; do
        # Collect all bins for this binner across all samples
        mkdir -p combined_\${binner}
        
        total_bins=0
        for meta_id in ${meta_ids}; do
            if [ -d "batch_bins/\${meta_id}_\${binner}" ]; then
                bin_count=\$(find -L batch_bins/\${meta_id}_\${binner} -maxdepth 1 -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) 2>/dev/null | wc -l)
                if [ "\$bin_count" -gt 0 ]; then
                    cp batch_bins/\${meta_id}_\${binner}/* combined_\${binner}/ 2>/dev/null || true
                    total_bins=\$((total_bins + bin_count))
                fi
            fi
        done
        
        if [ "\$total_bins" -eq 0 ]; then
            echo "WARNING: No genome bins found for \$binner across all samples. Creating empty reports."
            
            # Create empty reports for each sample
            for meta_id in ${meta_ids}; do
                echo -e "Name\\tCompleteness\\tContamination\\tCompleteness_Model_Used\\tTranslation_Table_Used\\tCoding_Density\\tContig_N50\\tAverage_Gene_Length\\tGenome_Size\\tGC_Content\\tTotal_Coding_Sequences\\tAdditional_Notes" > \${meta_id}_\${binner}_quality_report.tsv
            done
        else
            echo "Found \$total_bins genome bins for \$binner across all samples. Running CheckM2 assessment..."
            
            checkm2 predict \\
                --threads ${task.cpus} \\
                --input combined_\${binner} \\
                --output-directory checkm2_\${binner}_output \\
                --database_path ${checkm_db} \\
                -x .fa
            
            # Split the combined report by sample prefix
            if [ -f checkm2_\${binner}_output/quality_report.tsv ]; then
                # Process the combined report and split by sample
                python3 <<-PYEOF
				import sys
				import os
				
				binner = "\${binner}"
				meta_ids = "${meta_ids}".split()
				
				# Read the combined report
				with open(f"checkm2_{binner}_output/quality_report.tsv", 'r') as f:
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
				    with open(f"{meta_id}_{binner}_quality_report.tsv", 'w') as f:
				        f.write(header)
				        if sample_reports[meta_id]:
				            f.writelines(sample_reports[meta_id])
				PYEOF
            fi
        fi
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version 2>&1 | sed 's/checkm2: version //')
    END_VERSIONS
    """
    
    stub:
    def meta_ids = meta_list.collect { m -> m instanceof Map ? (m.id ?: m['id']) : m.toString() }.join(' ')
    """
    # Create stub quality report files for each sample and binner
    for meta_id in ${meta_ids}; do
        for binner in metabat semibin comebin metawrap; do
            cat > \${meta_id}_\${binner}_quality_report.tsv << 'EOF'
Name	Completeness	Contamination
bin.1	95.5	2.1
bin.2	87.3	1.5
EOF
        done
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: 1.1.0
    END_VERSIONS
    """
}
