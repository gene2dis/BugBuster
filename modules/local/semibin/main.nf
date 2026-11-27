/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SEMIBIN Binning Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenomic binning using SemiBin2
    
    Supports both per-sample and co-assembly modes:
    - Pre-trained model (single_easy_bin --environment): every per-sample
      assembly (params.semibin_env_model) and, by default, the co-assembly
      (params.semibin_coassembly_env_model) on its pooled BAM
    - Self-trained model (generate features → train_self → bin_short): the
      co-assembly with params.semibin_coassembly_mode = 'train_self' (Phase
      7a: 24.8 CPU-hours against 0.65 for the pre-trained model; after
      Binette both gave 36 bins, 25 and 23 of them high quality)
    
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
    
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/semibin:2.6.0--pyhdfd78af_0' :
        'quay.io/biocontainers/semibin:2.6.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs), path(bam)

    output:
    tuple val(meta), path("${meta.id}_semibin_output_bins"), emit: bins
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Extra tool options from a -c config (withName ... ext.args), as in every
    // other tool module (review R2-14); none set by default
    def args = task.ext.args ?: ''
    prefix = meta.id
    def is_coassembly = meta.coassembly as boolean
    def train_self = is_coassembly && params.semibin_coassembly_mode == 'train_self'
    // --semibin_env_model is null unless the user sets it (main.nf warns): human_gut
    def env_model = is_coassembly ? params.semibin_coassembly_env_model : (params.semibin_env_model ?: 'human_gut')
    // A fixed seed: unseeded SemiBin2 gives a bin more or less between reruns
    def seed = params.semibin_random_seed
    
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
    
    if (train_self) {
        // Co-assembly, self-trained model: multi-step workflow
        """
        set -euo pipefail
        
        # No contig of 1,000 bp at all: skip SemiBin2 without running it. A
        # contig set with some >= 1,000 bp but none >= 4,000 bp runs and is
        # caught from SemiBin2's log (no must-link pairs) by
        # bin/semibin_collect_bins.sh, which writes SKIPPED.txt when SemiBin2
        # found no bins and keeps any other failure's exit status
        min_length=1000
        long_contigs=\$(awk -v min=\$min_length '/^>/ {if (seqlen >= min) count++; seqlen=0; next} {seqlen += length(\$0)} END {if (seqlen >= min) count++; print count+0}' ${contigs})
        
        if [ "\$long_contigs" -eq 0 ]; then
            echo "WARNING: All contigs < \${min_length}bp. Skipping SemiBin for ${prefix}." >&2
            mkdir -p ${prefix}_semibin_output_bins
            echo "All contigs < \${min_length}bp - SemiBin skipped" > ${prefix}_semibin_output_bins/SKIPPED.txt
        else
            # Self-trained co-assembly workflow. Each step's output goes to
            # semibin.log and a failure skips the later steps: SemiBin2 stops
            # already at the features step when no contig reaches 2,500 or
            # 4,000 bp, which semibin_collect_bins.sh then reports as no bins
            status=0
            SemiBin2 generate_sequence_features_single \\
                     -i ${contigs} \\
                     -b ${bam} \\
                     -o contig_output \\
                     --threads ${task.cpus} \\
                     > semibin.log 2>&1 || status=\$?

            if [ "\$status" -eq 0 ]; then
                SemiBin2 train_self \\
                         --data contig_output/data.csv \\
                         --data-split contig_output/data_split.csv \\
                         -o contig_output \\
                         --random-seed ${seed} \\
                         --threads ${task.cpus} \\
                         >> semibin.log 2>&1 || status=\$?
            fi

            if [ "\$status" -eq 0 ]; then
                # SemiBin2 2.6.0 saves the model as model.pt (2.1.0: model.h5)
                model=contig_output/model.pt
                [ -f "\$model" ] || model=contig_output/model.h5
                SemiBin2 bin_short \\
                         ${args} \\
                         -i ${contigs} \\
                         --model "\$model" \\
                         --random-seed ${seed} \\
                         --data contig_output/data.csv \\
                         -o output \\
                         --compression none \\
                         --threads ${task.cpus} \\
                         >> semibin.log 2>&1 || status=\$?
            fi
            cat semibin.log >&2
            semibin_collect_bins.sh "\$status" output/output_bins ${prefix}_semibin_output_bins semibin.log ${prefix}
            rm -rf output contig_output semibin.log
        fi
        
        chmod 777 -R ${prefix}_semibin_output_bins
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            semibin: \$(SemiBin2 --version 2>&1 | sed 's/SemiBin2 //')
        END_VERSIONS
        """
    } else {
        // Pre-trained model: single easy bin
        """
        set -euo pipefail
        
        # No contig of 1,000 bp at all: skip SemiBin2 without running it. A
        # contig set with some >= 1,000 bp but none >= 4,000 bp runs and is
        # caught from SemiBin2's log (no must-link pairs) by
        # bin/semibin_collect_bins.sh, which writes SKIPPED.txt when SemiBin2
        # found no bins and keeps any other failure's exit status
        min_length=1000
        long_contigs=\$(awk -v min=\$min_length '/^>/ {if (seqlen >= min) count++; seqlen=0; next} {seqlen += length(\$0)} END {if (seqlen >= min) count++; print count+0}' ${contigs})
        
        if [ "\$long_contigs" -eq 0 ]; then
            echo "WARNING: All contigs < \${min_length}bp. Skipping SemiBin for ${prefix}." >&2
            mkdir -p ${prefix}_semibin_output_bins
            echo "All contigs < \${min_length}bp - SemiBin skipped" > ${prefix}_semibin_output_bins/SKIPPED.txt
        else
            status=0
            SemiBin2 single_easy_bin \\
                     ${args} \\
                     -i ${contigs} \\
                     -b ${bam} \\
                     -o ${prefix}_semibin_bins \\
                     --environment ${env_model} \\
                     --random-seed ${seed} \\
                     --compression none \\
                     --threads ${task.cpus} \\
                     > semibin.log 2>&1 || status=\$?
            cat semibin.log >&2
            semibin_collect_bins.sh "\$status" ${prefix}_semibin_bins/output_bins ${prefix}_semibin_output_bins semibin.log ${prefix}
            rm -rf ${prefix}_semibin_bins semibin.log
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
        semibin: 2.6.0
    END_VERSIONS
    """
}
