#!/bin/bash
#
# Generate the committed BAM fixtures (written for audit item #22, the old
# per-bin BEDTOOLS awk division by zero for bins with no aligned reads; now
# read by tests/modules/bedtools_bin_depth.nf.test and
# tests/modules/featurecounts_genes.nf.test).
#
# Creates, under tests/data/bam/:
#   - test_bin.1_all_reads.bam  two 50 bp reads aligned to a 1 kb contig
#   - test_bin.2_all_reads.bam  @SQ header but no aligned reads
#     (genomeCoverageBed still emits zero-depth rows for this shape)
#   - test_bin.3_all_reads.bam  no reference sequences at all (empty bin FASTA
#     upstream) -> genomeCoverageBed emits nothing, which crashed the old awk
#     chain with "Division by zero"
#   - test_bin.4_all_reads.bam  two contigs: 1 kb with two reads, 100 bp with
#     none (per-contig mean 0.05); genomeCoverageBed's "genome" summary rows
#     must not count as a third contig or enter the mean
#
# Filenames follow the per-bin naming of the removed BOWTIE2_SAMTOOLS_DEPTH
# ({sample}_bin.N_all_reads.bam); no test parses them any more.
#
# Requirements: docker.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUT_DIR="${SCRIPT_DIR}/../data/bam"
mkdir -p "${OUT_DIR}"
OUT_DIR="$( cd "${OUT_DIR}" && pwd )"

SAMTOOLS_IMG="quay.io/biocontainers/samtools:1.17--h00cdaf9_0"

READ_SEQ="$(printf 'ACGT%.0s' {1..12})AC"   # 50 bp
READ_QUAL="$(printf 'I%.0s' {1..50})"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

# Bin 1: two aligned reads on a 1 kb contig
cat > "${TMP_DIR}/bin1.sam" <<EOF
@HD	VN:1.6	SO:coordinate
@SQ	SN:contig_1	LN:1000
read_1	0	contig_1	1	60	50M	*	0	0	${READ_SEQ}	${READ_QUAL}
read_2	0	contig_1	101	60	50M	*	0	0	${READ_SEQ}	${READ_QUAL}
EOF

# Bin 2: @SQ header, no aligned reads
cat > "${TMP_DIR}/bin2.sam" <<EOF
@HD	VN:1.6	SO:coordinate
@SQ	SN:contig_2	LN:1000
EOF

# Bin 3: no reference sequences (the division-by-zero shape)
cat > "${TMP_DIR}/bin3.sam" <<EOF
@HD	VN:1.6	SO:coordinate
EOF

# Bin 4: two contigs of different length, reads on the first only
cat > "${TMP_DIR}/bin4.sam" <<EOF
@HD	VN:1.6	SO:coordinate
@SQ	SN:contig_a	LN:1000
@SQ	SN:contig_b	LN:100
read_1	0	contig_a	1	60	50M	*	0	0	${READ_SEQ}	${READ_QUAL}
read_2	0	contig_a	101	60	50M	*	0	0	${READ_SEQ}	${READ_QUAL}
EOF

docker run --rm -u "$(id -u):$(id -g)" \
    -v "${TMP_DIR}:/work" -w /work \
    "${SAMTOOLS_IMG}" sh -c '
        samtools view -b bin1.sam > test_bin.1_all_reads.bam
        samtools view -b bin2.sam > test_bin.2_all_reads.bam
        samtools view -b bin3.sam > test_bin.3_all_reads.bam
        samtools view -b bin4.sam > test_bin.4_all_reads.bam
    '

mv "${TMP_DIR}"/test_bin.*_all_reads.bam "${OUT_DIR}/"
echo "Fixtures written under ${OUT_DIR}"
