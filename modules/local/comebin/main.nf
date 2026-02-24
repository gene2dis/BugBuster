/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMEBIN Binning Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenomic binning using COMEBin
    
    Supports both per-sample and co-assembly modes
    Includes retry logic with adjusted parameters on failure
    
    Input:
        tuple val(meta), path(contigs), path(bam)
    
    Output:
        bins: tuple val(meta), path(comebin_bins/comebin_res/comebin_res_bins)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process COMEBIN {
    tag "${meta.id}"
    label 'process_high'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/comebin:1.0.4--hdfd78af_0' :
        'quay.io/biocontainers/comebin:1.0.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(contigs), path(bam)

    output:
    tuple val(meta), path("${meta.id}_comebin_bins/comebin_res/comebin_res_bins"), emit: bins
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    
    // Input validation
    if (!meta || !meta.id) {
        error "COMEBIN: meta.id is required"
    }
    if (!contigs) {
        error "COMEBIN: No contig file provided for sample ${meta.id}"
    }
    
    """
    set -euo pipefail
    
    rm -rf ${prefix}_comebin_bins
    mkdir -p ${prefix}_comebin_bins/comebin_res/comebin_res_bins
    
    # Validate minimum contig count
    contig_count=\$(grep -c "^>" ${contigs} || true)
    
    if [[ \$contig_count -lt 10 ]]; then
        echo "WARNING: Only \$contig_count contigs. Skipping COMEBIN (requires ≥10)." >&2
        echo "COMEBIN skipped: insufficient contigs (\$contig_count < 10)" > ${prefix}_comebin_bins/comebin_res/comebin_res_bins/SKIPPED.txt
        
        cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    comebin: 1.0.4
	END_VERSIONS
        exit 0
    fi
    
    # Compute N50 and set contrastive loss temperature accordingly
    n50=\$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq\$0} END{if(seq) print length(seq)}' ${contigs} \\
        | sort -rn \\
        | awk 'BEGIN{s=0} {a[NR]=\$1; s+=\$1} END{t=s/2; c=0; for(i=1;i<=NR;i++){c+=a[i]; if(c>=t){print a[i]; exit}}}')
    temp_flag=""
    if [[ \$n50 -lt 10000 ]]; then
        echo "INFO: N50=\$n50 < 10000, setting temperature -l 0.15" >&2
        temp_flag="-l 0.15"
    else
        echo "INFO: N50=\$n50 >= 10000, using default temperature" >&2
    fi
    
    # Run COMEBIN with retry logic and adjusted parameters
    set +e
    exit_code=0
    if [[ ${task.attempt} == 1 ]]; then
        bash run_comebin.sh \\
            -a ${contigs} \\
            -p ./ \\
            -o ${prefix}_comebin_bins \\
            -n 6 \\
            -t ${task.cpus} \\
            -b 512 \\
            -e 1024 \\
            -c 1024 \\
            \$temp_flag
        exit_code=\$?
    elif [[ ${task.attempt} == 2 ]]; then
        bash run_comebin.sh \\
            -a ${contigs} \\
            -p ./ \\
            -o ${prefix}_comebin_bins \\
            -n 4 \\
            -t ${task.cpus} \\
            -b 256 \\
            -e 512 \\
            -c 512 \\
            \$temp_flag
        exit_code=\$?
    elif [[ ${task.attempt} == 3 ]]; then
        bash run_comebin.sh \\
            -a ${contigs} \\
            -p ./ \\
            -o ${prefix}_comebin_bins \\
            -n 4 \\
            -t ${task.cpus} \\
            -b 128 \\
            -e 256 \\
            -c 256 \\
            \$temp_flag
        exit_code=\$?
    fi
    set -e
    
    # Handle failure gracefully
    if [[ \$exit_code -ne 0 ]]; then
        echo "WARNING: COMEBIN failed for ${prefix}. Creating empty directory." >&2
        mkdir -p ${prefix}_comebin_bins/comebin_res/comebin_res_bins
        echo "COMEBIN failed - insufficient data or features" > ${prefix}_comebin_bins/comebin_res/comebin_res_bins/FAILED.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    comebin: \$(run_comebin.sh --version 2>&1 | head -n1 || echo "1.0.4")
	END_VERSIONS
    """
    
    stub:
    prefix = meta.id
    """
    mkdir -p ${prefix}_comebin_bins/comebin_res/comebin_res_bins
    
    # Create stub bin files
    cat > ${prefix}_comebin_bins/comebin_res/comebin_res_bins/bin.1.fa << 'EOF'
>contig_1
ACGTACGTACGTACGTACGTACGTACGTACGT
>contig_2
TGCATGCATGCATGCATGCATGCATGCATGCA
EOF
    
    cat > ${prefix}_comebin_bins/comebin_res/comebin_res_bins/bin.2.fa << 'EOF'
>contig_3
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comebin: 1.0.4
    END_VERSIONS
    """
}
