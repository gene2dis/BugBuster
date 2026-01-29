/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SEMIBIN Binning Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenomic binning using SemiBin2
    
    Supports both per-sample and co-assembly modes:
    - Per-sample: Uses single_easy_bin workflow
    - Co-assembly: Uses multi-step workflow (generate features → train → bin)
    
    Input:
        tuple val(meta), path(contigs), path(bam)
    
    Output:
        bins: tuple val(meta), path(semibin_output_bins)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process SEMIBIN {
    tag "${meta.id}"
    label 'process_medium'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/semibin:2.1.0--pyhdfd78af_0' :
        'quay.io/biocontainers/semibin:2.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs), path(bam)

    output:
    tuple val(meta), path("${meta.id}_semibin_output_bins"), emit: bins
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    def env_model = params.semibin_env_model
    def is_coassembly = (meta.id == 'coassembly')
    
    // Input validation
    if (!meta || !meta.id) {
        error "SEMIBIN: meta.id is required"
    }
    if (!contigs) {
        error "SEMIBIN: No contig file provided for sample ${meta.id}"
    }
    if (!bam) {
        error "SEMIBIN: No BAM file(s) provided for sample ${meta.id}"
    }
    
    if (is_coassembly) {
        // Co-assembly mode: multi-step workflow
        """
        set -euo pipefail
        
        # Validate minimum contig length
        min_length=1000
        long_contigs=\$(awk -v min=\$min_length '/^>/ {if (seqlen >= min) count++; seqlen=0; next} {seqlen += length(\$0)} END {if (seqlen >= min) count++; print count+0}' ${contigs})
        
        if [ "\$long_contigs" -eq 0 ]; then
            echo "WARNING: All contigs < \${min_length}bp. Skipping SemiBin for ${prefix}." >&2
            mkdir -p ${prefix}_semibin_output_bins
            echo "All contigs < \${min_length}bp - SemiBin skipped" > ${prefix}_semibin_output_bins/SKIPPED.txt
        else
            # Multi-sample co-assembly workflow
            SemiBin2 generate_sequence_features_single \\
                     -i ${contigs} \\
                     -b ${bam} \\
                     -o contig_output \\
                     --threads ${task.cpus}
            
            SemiBin2 train_self \\
                     --data contig_output/data.csv \\
                     --data-split contig_output/data_split.csv \\
                     -o contig_output \\
                     --threads ${task.cpus}
            
            SemiBin2 bin_short \\
                     -i ${contigs} \\
                     --model contig_output/model.h5 \\
                     --data contig_output/data.csv \\
                     -o output \\
                     --compression none \\
                     --threads ${task.cpus}
            
            mv output/output_bins ${prefix}_semibin_output_bins
            rm -rf output contig_output
        fi
        
        chmod 777 -R ${prefix}_semibin_output_bins
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            semibin: \$(SemiBin2 --version 2>&1 | sed 's/SemiBin2 //')
        END_VERSIONS
        """
    } else {
        // Per-sample mode: single easy bin
        """
        set -euo pipefail
        
        # Validate minimum contig length
        min_length=1000
        long_contigs=\$(awk -v min=\$min_length '/^>/ {if (seqlen >= min) count++; seqlen=0; next} {seqlen += length(\$0)} END {if (seqlen >= min) count++; print count+0}' ${contigs})
        
        if [ "\$long_contigs" -eq 0 ]; then
            echo "WARNING: All contigs < \${min_length}bp. Skipping SemiBin for ${prefix}." >&2
            mkdir -p ${prefix}_semibin_output_bins
            echo "All contigs < \${min_length}bp - SemiBin skipped" > ${prefix}_semibin_output_bins/SKIPPED.txt
        else
            SemiBin2 single_easy_bin \\
                     -i ${contigs} \\
                     -b ${bam} \\
                     -o ${prefix}_semibin_bins \\
                     --environment ${env_model} \\
                     --compression none \\
                     --threads ${task.cpus}
            
            mv ${prefix}_semibin_bins/output_bins ${prefix}_semibin_output_bins
            rm -rf ${prefix}_semibin_bins
        fi
        
        chmod 777 -R ${prefix}_semibin_output_bins
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            semibin: \$(SemiBin2 --version 2>&1 | sed 's/SemiBin2 //')
        END_VERSIONS
        """
    }
    
    stub:
    prefix = meta.id
    """
    mkdir -p ${prefix}_semibin_output_bins
    
    # Create stub bin files
    cat > ${prefix}_semibin_output_bins/bin.1.fa << 'EOF'
>contig_1
ACGTACGTACGTACGTACGTACGTACGTACGT
>contig_2
TGCATGCATGCATGCATGCATGCATGCATGCA
EOF
    
    cat > ${prefix}_semibin_output_bins/bin.2.fa << 'EOF'
>contig_3
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        semibin: 2.1.0
    END_VERSIONS
    """
}
