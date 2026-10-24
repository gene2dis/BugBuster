/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BBNORM Read Normalization Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Digital normalization of the assembler input with BBNorm
    (--assembly_read_reduction normalize). Only MEGAHIT gets these reads:
    mapping, coverage, binning and gene counting use the full clean reads.

    Supports both per-sample and co-assembly modes (the pooled reads of every
    sample, meta.id 'coassembly', are normalized together in one call)

    Input:
        tuple val(meta), path(reads)

    Output:
        reads: tuple val(meta), path(normalized reads: R1, R2(, singleton))
        report: path(read_reduction.tsv)
        log: path(bbnorm.log)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process BBNORM {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:40.05--h09cc210_0' :
        'quay.io/biocontainers/bbmap:40.05--h09cc210_0' }"

    // Publishing is configured in config/modules.config (withName: 'BBNORM')

    input:
    tuple val(meta), path(reads)

    output:
    // Sorted by name: R1, R2, then the singletons, the order MEGAHIT's
    // per-sample branch reads them in
    tuple val(meta), path("${meta.id}_norm_*.fastq*")  , emit: reads
    path "${meta.id}_read_reduction.tsv"               , emit: report
    path "${meta.id}_bbnorm.log"                       , emit: log
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    def args = task.ext.args ?: ''

    // Input validation
    if (!meta || !meta.id) {
        error "BBNORM: meta.id is required"
    }
    if (!reads || reads.size() < 2) {
        error "BBNORM: paired reads (R1 and R2) are required for sample ${meta.id}"
    }

    // bbnorm.sh sizes its heap from the host's free memory, not the task's
    // limit, and fills whatever heap it gets (its k-mer tables grow to it):
    // 85 % of the task memory leaves room for the JVM outside the heap
    def heap = task.memory ? "-Xmx${(task.memory.toMega() * 0.85) as long}m" : ''

    // Co-assembly pools every sample's reads; bin/read_role.sh sorts them by
    // file name, as MEGAHIT and the aligners do. A single sample's reads
    // arrive positionally (raw input names when QC is off).
    def classify_reads = meta.id == 'coassembly'
        ? """for read_file in ${reads.join(' ')}; do
        role=\$(read_role.sh "\$read_file") || exit 1
        case "\$role" in
            R1) R1_files+=("\$read_file") ;;
            R2) R2_files+=("\$read_file") ;;
            S)  S_files+=("\$read_file") ;;
        esac
    done"""
        : "R1_files=(${reads[0]}); R2_files=(${reads[1]})${reads.size() > 2 ? "; S_files=(${reads[2]})" : ''}"

    """
    set -euo pipefail

    R1_files=()
    R2_files=()
    S_files=()
    ${classify_reads}
    if [ \${#R1_files[@]} -eq 0 ] || [ \${#R1_files[@]} -ne \${#R2_files[@]} ]; then
        echo "ERROR: BBNORM needs matching R1 and R2 files; got: ${reads}" >&2
        exit 1
    fi

    # Pairs in and kept, counted from the files (BBNorm logs reads per input
    # file pair and per pass). The input is counted while BBNorm runs.
    count_reads() {
        case "\$1" in
            *.gz) gzip -cd "\$1" ;;
            *)    cat "\$1" ;;
        esac | awk 'END { print NR / 4 }'
    }
    count_pairs_in() {
        local n=0 c f
        for f in "\${R1_files[@]}"; do
            c=\$(count_reads "\$f")
            n=\$(( n + c ))
        done
        echo "\$n" > .pairs_in
    }
    count_pairs_in &
    count_pid=\$!

    # BBNorm reads its input more than once, so no pipes: several files go in
    # as comma-separated lists, read one pair of files after the other
    in1=\$(IFS=,; echo "\${R1_files[*]}")
    in2=\$(IFS=,; echo "\${R2_files[*]}")

    # Its log is published; on failure its end goes to stderr for the report
    if ! bbnorm.sh ${heap} \\
        in1="\$in1" \\
        in2="\$in2" \\
        out1=${prefix}_norm_R1.fastq.gz \\
        out2=${prefix}_norm_R2.fastq.gz \\
        threads=${task.cpus} \\
        zl=4 \\
        ${args} \\
        2> ${prefix}_bbnorm.log; then
        kill \$count_pid 2> /dev/null || true
        tail -n 50 ${prefix}_bbnorm.log >&2
        exit 1
    fi

    # Singletons go to the assembler unchanged (links, compression kept):
    # BBNorm cannot write paired and unpaired reads in one call, and they are
    # a small share of the reads. Several (co-assembly) stay separate files.
    singleton_suffix() {
        case "\$1" in
            *.gz)  echo fastq.gz ;;
            *.bz2) echo fastq.bz2 ;;
            *)     echo fastq ;;
        esac
    }
    if [ \${#S_files[@]} -eq 1 ]; then
        ln -s "\${S_files[0]}" ${prefix}_norm_singleton.\$(singleton_suffix "\${S_files[0]}")
    else
        for i in "\${!S_files[@]}"; do
            ln -s "\${S_files[i]}" ${prefix}_norm_singleton_\$(( i + 1 )).\$(singleton_suffix "\${S_files[i]}")
        done
    fi

    wait \$count_pid
    pairs_in=\$(cat .pairs_in)
    pairs_kept=\$(count_reads ${prefix}_norm_R1.fastq.gz)

    printf 'sample_id\\tmethod\\tpairs_in\\tpairs_kept\\tpercent_kept\\n' > ${prefix}_read_reduction.tsv
    awk -v s=${prefix} -v i="\$pairs_in" -v k="\$pairs_kept" \\
        'BEGIN { printf "%s\\tbbnorm\\t%d\\t%d\\t%.2f\\n", s, i, k, (i > 0 ? 100 * k / i : 0) }' \\
        >> ${prefix}_read_reduction.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh 2>&1 | head -n 1)
    END_VERSIONS
    """

    stub:
    prefix = meta.id
    def classify_reads = meta.id == 'coassembly'
        ? """for read_file in ${reads.join(' ')}; do
        if [ "\$(read_role.sh "\$read_file")" = S ]; then S_files+=("\$read_file"); fi
    done"""
        : (reads.size() > 2 ? "S_files=(${reads[2]})" : '')
    """
    echo "@stub/1" | gzip > ${prefix}_norm_R1.fastq.gz
    echo "@stub/2" | gzip > ${prefix}_norm_R2.fastq.gz
    S_files=()
    ${classify_reads}
    if [ \${#S_files[@]} -gt 0 ]; then
        echo "@stub" | gzip > ${prefix}_norm_singleton.fastq.gz
    fi
    printf 'sample_id\\tmethod\\tpairs_in\\tpairs_kept\\tpercent_kept\\n${prefix}\\tbbnorm\\t10\\t8\\t80.00\\n' > ${prefix}_read_reduction.tsv
    touch ${prefix}_bbnorm.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: 40.05
    END_VERSIONS
    """
}
