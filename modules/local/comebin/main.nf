/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMEBIN Binning Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenomic binning using COMEBin

    Supports both per-sample and co-assembly modes. Runs on CPU only (the
    image's PyTorch has no CUDA). Each retry uses smaller batch settings; a
    failure on the last attempt leaves FAILED.txt instead of failing the run

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
        'https://depot.galaxyproject.org/singularity/comebin:1.1.0--hdfd78af_2' :
        'quay.io/biocontainers/comebin:1.1.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(contigs), path(bam)

    output:
    tuple val(meta), path("${meta.id}_comebin_bins/comebin_res/comebin_res_bins"), emit: bins
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    // Batch settings per attempt: -n views, -b batch size, -e/-c embedding sizes
    def settings = task.attempt == 1 ? '-n 6 -b 512 -e 1024 -c 1024' :
        task.attempt == 2 ? '-n 4 -b 256 -e 512 -c 512' :
        '-n 4 -b 128 -e 256 -c 256'
    // A failed attempt fails the task so the retry runs with smaller batches;
    // only the last one writes FAILED.txt and lets the rest of the run go on
    def last_attempt = task.attempt > (task.maxRetries ?: 0)

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
    bins_dir=${prefix}_comebin_bins/comebin_res/comebin_res_bins

    write_versions() {
        printf '"%s":\\n    comebin: %s\\n' "${task.process}" "\$(run_comebin.sh -V 2>&1 | sed 's/.*: //')" > versions.yml
    }

    # Validate minimum contig count (awk always exits 0 and always prints a
    # number; 'grep -c ... || true' left the variable empty when grep failed)
    contig_count=\$(awk '/^>/{n++} END{print n+0}' ${contigs})

    if [[ \$contig_count -lt 10 ]]; then
        echo "WARNING: Only \$contig_count contigs. Skipping COMEBIN (requires ≥10)." >&2
        echo "COMEBIN skipped: insufficient contigs (\$contig_count < 10)" > "\$bins_dir/SKIPPED.txt"
        write_versions
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

    # run_comebin.sh calls 'realpath -e --' and 'realpath -m --', which the
    # image's busybox realpath does not support: a Python realpath written
    # into the task directory comes first on PATH (nothing is written into
    # the image)
    mkdir -p realpath_shim
    cat > realpath_shim/realpath <<'PYEOF'
    #!/usr/bin/env python3
    import os
    import sys
    args = sys.argv[1:]
    must_exist = '-e' in args
    status = 0
    for arg in [a for a in args if a not in ('-e', '-m', '--')]:
        path = os.path.realpath(arg)
        if must_exist and not os.path.exists(path):
            sys.stderr.write('realpath: ' + arg + ': No such file or directory\\n')
            status = 1
            continue
        print(path)
    sys.exit(status)
    PYEOF
    chmod +x realpath_shim/realpath
    export PATH="\$PWD/realpath_shim:\$PATH"

    # -d cpu: the default device is cuda, and the image's PyTorch has no CUDA.
    # -s: a fixed seed, so reruns give the same bins
    exit_code=0
    bash run_comebin.sh \\
        -a ${contigs} \\
        -p ./ \\
        -o ${prefix}_comebin_bins \\
        -t ${task.cpus} \\
        -d cpu \\
        -s ${params.comebin_random_seed} \\
        ${settings} \\
        \$temp_flag || exit_code=\$?

    if [[ \$exit_code -ne 0 ]]; then
        if [[ "${last_attempt}" != "true" ]]; then
            echo "ERROR: COMEBIN failed for ${prefix} (exit \$exit_code) on attempt ${task.attempt}; retrying with smaller batches." >&2
            exit \$exit_code
        fi
        echo "WARNING: COMEBIN failed for ${prefix} on its last attempt. Creating empty directory." >&2
        rm -rf ${prefix}_comebin_bins
        mkdir -p "\$bins_dir"
        echo "COMEBIN failed on every attempt (last exit status \$exit_code)" > "\$bins_dir/FAILED.txt"
    else
        # COMEBin 1.1.0 keeps the whole contig header (e.g. MEGAHIT's
        # 'flag=... multi=... len=...') in its bins; refiners and the bin
        # depth step match contigs by name, so keep the name only
        find "\$bins_dir" -name '*.fa' -exec sed -i '/^>/s/[[:space:]].*//' {} +
    fi

    write_versions
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
        comebin: 1.1.0
    END_VERSIONS
    """
}
