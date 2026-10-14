/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BWAMEM3_SAMTOOLS Alignment Module (--read_aligner bwamem3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Align reads to assembled contigs (or to all refined bins) with bwa-mem3
    and write one coordinate-sorted BAM per sample.

    The index is its own task (BWAMEM3_INDEX), so one co-assembly index serves
    the pooled alignment and every per-sample alignment instead of being
    rebuilt in each of them. BINS_REFERENCE concatenates the refined bins into
    one reference plus a contig -> bin table, so co-assembly bin depth maps
    each sample once: bwa mem aligns locally, and mapping against one bin at
    a time would let reads from other bins pile up as clipped partial hits.

    Alignments pass bin/filter_alignments.sh before sorting: a primary
    alignment needs >= ext.min_aligned_length aligned bases and an identity
    >= ext.min_identity, otherwise the read is kept as unmapped (samtools
    fixmate then updates its mate). bwa mem's local hits from genomes the
    reference does not contain otherwise inflate gene counts and bin depth
    in poorly covered samples (docs/benchmarks/bwamem3_vs_bowtie2.md).

    Input:
        BWAMEM3_INDEX:    tuple val(meta), path(fasta)
        BWAMEM3_SAMTOOLS: tuple val(meta), path(reads), path(contigs), path(index)
        BINS_REFERENCE:   path(bins)

    Output:
        BWAMEM3_INDEX:    index: tuple val(meta), path(index dir)
        BWAMEM3_SAMTOOLS: contigs_and_bam, bam_only, empty_contig_report
        BINS_REFERENCE:   reference: tuple path(all_bins.fa), path(contig_to_bin.tsv)
----------------------------------------------------------------------------------------
*/

process BWAMEM3_INDEX {
    tag "${meta.id}"
    // No resource label with a memory setting: withLabel would override the
    // size-based memory below. process_index only sets cpus and time.
    label 'process_index'
    container 'community.wave.seqera.io/library/bwa-mem3_samtools:65dfdfc4140db421'

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
        bwamem3: 0.14.0
    END_VERSIONS
    """
}

process BWAMEM3_SAMTOOLS {
    tag "${meta.id}"
    label 'process_high'
    container 'community.wave.seqera.io/library/bwa-mem3_samtools:65dfdfc4140db421'

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

    // Co-assembly pools every sample's reads (sorted by sample, each sample
    // R1, R2[, singleton]); bin/read_role.sh sorts them by file name, as in
    // MEGAHIT and BOWTIE2_SAMTOOLS, so all three agree on the pairs. bwa mem
    // takes one file per mate, not bowtie2's comma lists, so each library is
    // mapped on its own and the BAMs are merged; that also keeps bwa's
    // insert-size estimate per library. A single sample's reads arrive
    // positionally (raw input names when QC is off).
    def classify_reads = meta.id == 'coassembly'
        ? """for f in ${reads.join(' ')}; do
        case "\$(read_role.sh "\$f")" in
            R1) R1+=("\$f") ;;
            R2) R2+=("\$f") ;;
            S)  S+=("\$f") ;;
        esac
    done"""
        : "R1=(${reads[0]}); R2=(${reads[1]})${reads.size() > 2 ? "; S=(${reads[2]})" : ''}"
    def stage = task.process.tokenize(':')[-1]

    """
    set -euo pipefail

    R1=(); R2=(); S=()
    ${classify_reads}
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
        bwamem3: 0.14.0
        samtools: 1.22.1
    END_VERSIONS
    """
}

process BINS_REFERENCE {
    tag "bins"
    label 'process_single'
    container 'community.wave.seqera.io/library/bwa-mem3_samtools:65dfdfc4140db421'

    input:
    path(bins)

    output:
    tuple path("all_bins.fa"), path("contig_to_bin.tsv"), emit: reference
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    # contig_to_bin.tsv: <contig> TAB <bin>, bins in glob order. A bin with no
    # sequences gets one row with an empty contig field, so it still reports
    # a 0-coverage row downstream. Bin ids are the FASTA basenames. A contig
    # already placed in an earlier bin is skipped with a warning, so the
    # reference has unique names (MetaWRAP-refined bins should never share
    # contigs; the old per-bin mapping tolerated it silently).
    : > all_bins.fa
    : > contig_to_bin.tsv
    for bin in ${bins}/*.fa; do
        [ -e "\${bin}" ] || continue
        bin_name=\$(basename "\${bin}" .fa)
        if ! grep -q '^>' "\${bin}"; then
            printf '\\t%s\\n' "\${bin_name}" >> contig_to_bin.tsv
            continue
        fi
        awk -v bin="\${bin_name}" '
            BEGIN {
                while ((getline line < "contig_to_bin.tsv") > 0) {
                    split(line, f, "\\t")
                    if (f[1] != "") seen[f[1]] = 1
                }
            }
            /^>/ {
                name = substr(\$1, 2)
                keep = !(name in seen)
                if (keep) {
                    seen[name] = 1
                    print name "\\t" bin > "contig_to_bin.part"
                } else {
                    print "WARNING: contig " name " of bin " bin " is already in another bin - skipped" > "/dev/stderr"
                }
            }
            keep { print }
        ' "\${bin}" >> all_bins.fa
        if [ -f contig_to_bin.part ]; then
            cat contig_to_bin.part >> contig_to_bin.tsv
            rm contig_to_bin.part
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk -W version 2>&1 | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    printf '>contig_1\\nACGT\\n' > all_bins.fa
    printf 'contig_1\\tbin.1\\n' > contig_to_bin.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: 1.3.4
    END_VERSIONS
    """
}
