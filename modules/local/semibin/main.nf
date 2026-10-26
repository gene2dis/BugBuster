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
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
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
    prefix = meta.id
    def is_coassembly = (meta.id == 'coassembly')
    def train_self = is_coassembly && params.semibin_coassembly_mode == 'train_self'
    def env_model = is_coassembly ? params.semibin_coassembly_env_model : params.semibin_env_model
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
        
        # SemiBin2 fails when it finds no bins: the marker-gene search of its
        # reclustering step gets an empty ORF file. An empty bin directory
        # with a note lets Binette and the rest of the sample go on; any
        # other failure keeps its exit status
        no_bins_or_move() {
            if [ "\$status" -eq 0 ]; then
                mv "\$1" ${prefix}_semibin_output_bins
            elif grep -q 'Number of bins prior to reclustering: 0' semibin.log; then
                echo "WARNING: SemiBin found no bins for ${prefix}." >&2
                mkdir -p ${prefix}_semibin_output_bins
                echo "SemiBin found no bins" > ${prefix}_semibin_output_bins/SKIPPED.txt
            else
                exit "\$status"
            fi
        }
        
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
                     --random-seed ${seed} \\
                     --threads ${task.cpus}
            
            # SemiBin2 2.6.0 saves the model as model.pt (2.1.0: model.h5)
            model=\$(ls contig_output/model.pt contig_output/model.h5 2>/dev/null | head -n 1)
            status=0
            SemiBin2 bin_short \\
                     -i ${contigs} \\
                     --model "\$model" \\
                     --random-seed ${seed} \\
                     --data contig_output/data.csv \\
                     -o output \\
                     --compression none \\
                     --threads ${task.cpus} \\
                     > semibin.log 2>&1 || status=\$?
            cat semibin.log >&2
            no_bins_or_move output/output_bins
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
        
        # SemiBin2 fails when it finds no bins: the marker-gene search of its
        # reclustering step gets an empty ORF file. An empty bin directory
        # with a note lets Binette and the rest of the sample go on; any
        # other failure keeps its exit status
        no_bins_or_move() {
            if [ "\$status" -eq 0 ]; then
                mv "\$1" ${prefix}_semibin_output_bins
            elif grep -q 'Number of bins prior to reclustering: 0' semibin.log; then
                echo "WARNING: SemiBin found no bins for ${prefix}." >&2
                mkdir -p ${prefix}_semibin_output_bins
                echo "SemiBin found no bins" > ${prefix}_semibin_output_bins/SKIPPED.txt
            else
                exit "\$status"
            fi
        }
        
        # Validate minimum contig length
        min_length=1000
        long_contigs=\$(awk -v min=\$min_length '/^>/ {if (seqlen >= min) count++; seqlen=0; next} {seqlen += length(\$0)} END {if (seqlen >= min) count++; print count+0}' ${contigs})
        
        if [ "\$long_contigs" -eq 0 ]; then
            echo "WARNING: All contigs < \${min_length}bp. Skipping SemiBin for ${prefix}." >&2
            mkdir -p ${prefix}_semibin_output_bins
            echo "All contigs < \${min_length}bp - SemiBin skipped" > ${prefix}_semibin_output_bins/SKIPPED.txt
        else
            status=0
            SemiBin2 single_easy_bin \\
                     -i ${contigs} \\
                     -b ${bam} \\
                     -o ${prefix}_semibin_bins \\
                     --environment ${env_model} \\
                     --random-seed ${seed} \\
                     --compression none \\
                     --threads ${task.cpus} \\
                     > semibin.log 2>&1 || status=\$?
            cat semibin.log >&2
            no_bins_or_move ${prefix}_semibin_bins/output_bins
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
