/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BWAMEM3_SAMTOOLS Alignment Module (--read_aligner bwamem3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Align reads to assembled contigs (a sample's own, or the co-assembly)
    with bwa-mem3 and write one coordinate-sorted BAM per sample.

    The index is its own task (BWAMEM3_INDEX), so one co-assembly index serves
    every per-sample alignment instead of being rebuilt in each of them.

    Alignments pass bin/filter_alignments.sh before sorting: a primary
    alignment needs >= ext.min_aligned_length aligned bases and an identity
    >= ext.min_identity, otherwise the read is kept as unmapped (samtools
    fixmate then updates its mate). bwa mem's local hits from genomes the
    reference does not contain otherwise inflate gene counts and depths
    in poorly covered samples (docs/benchmarks/bwamem3_vs_bowtie2.md).

    Input:
        BWAMEM3_INDEX:    tuple val(meta), path(fasta)
        BWAMEM3_SAMTOOLS: tuple val(meta), path(reads), path(contigs), path(index)

    Output:
        BWAMEM3_INDEX:    index: tuple val(meta), path(index dir)
        BWAMEM3_SAMTOOLS: contigs_and_bam, bam_only, empty_contig_report
----------------------------------------------------------------------------------------
*/

process BWAMEM3_INDEX {
    tag "${meta.id}"
    // No resource label with a memory setting: withLabel would override the
    // size-based memory below. process_index only sets cpus and time.
    label 'process_index'
    container 'community.wave.seqera.io/library/bwa-mem3_samtools:96d089551af9e9c3'

    // Peak ~8 bytes per base of doubled text below ~1.07 Gbp and ~12 above
    // (bwa-mem3 docs): 28 B/base covers both with headroom. The FASTA is
    // uncompressed, so its size is the base count (plus headers).
    memory { [2.GB, 280.MB * Math.ceil(fasta.size() / 10000000)].max() * task.attempt }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_bwamem3_index"), emit: index
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    """
    set -euo pipefail

    mkdir ${prefix}_bwamem3_index

    # bwa-mem3 cannot index a FASTA without sequences (empty contig set after
    # length filtering, no bins). Leave the index directory empty: the
    # mapping task sees the same empty FASTA and writes a header-only BAM.
    if grep -q '^>' ${fasta}; then
        bwa-mem3 index \\
            ${args} \\
            -t ${task.cpus} \\
            --max-memory ${task.memory.toMega()}M \\
            -p ${prefix}_bwamem3_index/${prefix} \\
            ${fasta}
    else
        echo "WARNING: ${fasta} has no sequences - leaving the index empty" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem3: \$(bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir ${prefix}_bwamem3_index
    touch ${prefix}_bwamem3_index/${prefix}.amb
    touch ${prefix}_bwamem3_index/${prefix}.ann
    touch ${prefix}_bwamem3_index/${prefix}.bwt.2bit.64
    touch ${prefix}_bwamem3_index/${prefix}.pac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem3: 0.15.0
    END_VERSIONS
    """
}

process BWAMEM3_SAMTOOLS {
    tag "${meta.id}"
    label 'process_high'
    container 'community.wave.seqera.io/library/bwa-mem3_samtools:96d089551af9e9c3'

    input:
    tuple val(meta), path(reads), path(contigs), path(index)

    output:
    tuple val(meta), path(contigs), path("${meta.id}_all_reads.bam"), emit: contigs_and_bam
    tuple val(meta), path("${meta.id}_all_reads.bam")              , emit: bam_only
    path "${meta.id}_empty_contig_report.txt", optional: true      , emit: empty_contig_report
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    // != null, not ?:, so that 0 (filter off) is honoured
    def min_identity = task.ext.min_identity != null ? task.ext.min_identity : 0.95
    def min_aligned_length = task.ext.min_aligned_length != null ? task.ext.min_aligned_length : 50
    def prefix = meta.id

    if (!reads || reads.size() < 2) {
        error "BWAMEM3_SAMTOOLS: paired reads are required for sample ${meta.id}"
    }

    // One sample per task, reads positional (raw input names when QC is
    // off). Under co-assembly each sample is mapped on its own against the
    // co-assembly and the pooled BAM is a samtools merge of these (issue
    // #28); bwa's insert-size estimate stays per library.
    def stage = task.process.tokenize(':')[-1]

    """
    set -euo pipefail

    R1=(${reads[0]}); R2=(${reads[1]}); S=(${reads.size() > 2 ? reads[2] : ''})
    if [ \${#R1[@]} -eq 0 ] || [ \${#R1[@]} -ne \${#R2[@]} ]; then
        echo "ERROR: could not pair R1/R2 files for ${prefix}: ${reads.join(' ')}" >&2
        exit 1
    fi

    # Contig set empty after filtering (BWAMEM3_INDEX left the index empty)
    if ! grep -q '^>' "${contigs}"; then
        echo "WARNING: Contig file is empty for ${prefix} - creating empty BAM" >&2

        printf '@HD\\tVN:1.6\\tSO:coordinate\\n' > empty.sam
        samtools view -b -o ${prefix}_all_reads.bam empty.sam

        printf 'sample_id\\t${prefix}\\n' > ${prefix}_empty_contig_report.txt
        printf 'stage\\t${stage}\\n' >> ${prefix}_empty_contig_report.txt
        printf 'empty_contig_detected\\tTRUE\\n' >> ${prefix}_empty_contig_report.txt
        printf 'timestamp\\t%s\\n' "\$(date -u +"%Y-%m-%d %H:%M:%S UTC")" >> ${prefix}_empty_contig_report.txt
        printf 'action_taken\\tCreated empty BAM file\\n' >> ${prefix}_empty_contig_report.txt
    else
        bwt=\$(find -L ${index} -maxdepth 1 -name '*.bwt.2bit.64' | sort | head -1)
        if [ -z "\${bwt}" ]; then
            echo "ERROR: no bwa-mem3 index (*.bwt.2bit.64) in ${index}" >&2
            exit 1
        fi
        idx=\${bwt%.bwt.2bit.64}

        # One coordinate-sorted BAM per library. filter_alignments.sh drops
        # secondary and supplementary records (downstream depth, binning and
        # featureCounts expect one record per read, as bowtie2 end-to-end
        # produced) and turns short or divergent hits into unmapped records;
        # fixmate updates their mates (bwa mem keeps a read's records
        # together, as fixmate needs). Unmapped records stay. -K fixes the
        # batch size so the output does not depend on task.cpus.
        map_reads() {
            local out=\$1
            shift
            bwa-mem3 mem ${args} -t ${task.cpus} -K 100000000 \${idx} "\$@" \\
                | filter_alignments.sh ${min_identity} ${min_aligned_length} \\
                | samtools fixmate -u - - \\
                | samtools sort ${args2} -@ ${task.cpus} -o "\${out}" -
        }

        bams=()
        for i in "\${!R1[@]}"; do
            bam="${prefix}_lib\${i}_paired.bam"
            map_reads "\${bam}" "\${R1[i]}" "\${R2[i]}"
            bams+=("\${bam}")
        done
        for i in "\${!S[@]}"; do
            bam="${prefix}_lib\${i}_singleton.bam"
            map_reads "\${bam}" "\${S[i]}"
            bams+=("\${bam}")
        done
        if [ \${#bams[@]} -eq 1 ]; then
            mv "\${bams[0]}" ${prefix}_all_reads.bam
        else
            samtools merge -@ ${task.cpus} -o ${prefix}_all_reads.bam "\${bams[@]}"
            rm -f "\${bams[@]}"
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem3: \$(bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    printf '@HD\\tVN:1.6\\tSO:coordinate\\n@SQ\\tSN:contig_1\\tLN:1000\\n' > header.sam
    samtools view -b -o ${prefix}_all_reads.bam header.sam 2>/dev/null || touch ${prefix}_all_reads.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem3: 0.15.0
        samtools: 1.24
    END_VERSIONS
    """
}
