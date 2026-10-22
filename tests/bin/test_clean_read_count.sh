#!/bin/bash
#
# Regression tests for bin/clean_read_count.sh, which BWAMEM3_DECONTAMINATE
# and BOWTIE2_DECONTAMINATE use to take the clean read count from the
# aligner's stderr summary instead of decompressing the output again
# (performance review E8). The logs below are real: samtools 1.24 and
# bowtie2 2.5.5 on 2 M read pairs of a real metagenome, where the counts
# matched `zcat | wc -l` of the written files, plus samtools on a
# constructed SAM with two orphan unmapped mates (6 processed, 2 singletons,
# 4 reads written to -1/-2), and bowtie2's whole summary for an empty input.
#
# Requirements: bash only.
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
COUNT="${REPO_DIR}/bin/clean_read_count.sh"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT
PASS=0
FAIL=0

expect() {
    local desc="$1" log="$2" mode="$3" want="$4" got
    got=$("${COUNT}" "${TMP_DIR}/${log}" "${mode}" 2>/dev/null) || got=""
    if [ "${got}" = "${want}" ]; then
        echo "✓ ${desc} -> ${want}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${desc}: expected '${want}', got '${got:-<error>}'"
        FAIL=$((FAIL + 1))
    fi
}

expect_error() {
    local desc="$1" log="$2" mode="$3"
    if "${COUNT}" "${TMP_DIR}/${log}" "${mode}" > /dev/null 2>&1; then
        echo "✗ ${desc}: expected an error, got exit 0"
        FAIL=$((FAIL + 1))
    else
        echo "✓ ${desc} -> error"
        PASS=$((PASS + 1))
    fi
}

cat > "${TMP_DIR}/samtools_paired.log" <<'EOF'
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 3756998 reads
EOF

cat > "${TMP_DIR}/samtools_orphans.log" <<'EOF'
[M::bam2fq_mainloop] discarded 2 singletons
[M::bam2fq_mainloop] processed 6 reads
EOF

cat > "${TMP_DIR}/samtools_single.log" <<'EOF'
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 298003 reads
EOF

cat > "${TMP_DIR}/bowtie2_paired.log" <<'EOF'
[WARNING] Failed to launch x86-64-v3 version, staying with default
[WARNING] Failed to launch x86-64-v3 version, staying with default
1892786 reads; of these:
  1892786 (100.00%) were paired; of these:
    1847541 (97.61%) aligned concordantly 0 times
    45245 (2.39%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
    ----
    1847541 pairs aligned concordantly 0 times; of these:
      1392 (0.08%) aligned discordantly 1 time
    ----
    1846149 pairs aligned 0 times concordantly or discordantly; of these:
      3692298 mates make up the pairs; of these:
        3656891 (99.04%) aligned 0 times
        35247 (0.95%) aligned exactly 1 time
        160 (0.00%) aligned >1 times
3.40% overall alignment rate
EOF

cat > "${TMP_DIR}/bowtie2_unpaired.log" <<'EOF'
[WARNING] Failed to launch x86-64-v3 version, staying with default
[WARNING] Failed to launch x86-64-v3 version, staying with default
300000 reads; of these:
  300000 (100.00%) were unpaired; of these:
    288005 (96.00%) aligned 0 times
    11995 (4.00%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
4.00% overall alignment rate
EOF

: > "${TMP_DIR}/empty.log"

# bowtie2 2.5.5 on an empty FASTQ (-U or -1/-2): no per-category lines
cat > "${TMP_DIR}/bowtie2_no_reads.log" <<'EOF'
[WARNING] Failed to launch x86-64-v3 version, staying with default
[WARNING] Failed to launch x86-64-v3 version, staying with default
0 reads
0.00% overall alignment rate
EOF

expect "samtools fastq, paired"               samtools_paired.log  samtools          3756998
expect "samtools fastq, orphan mates removed" samtools_orphans.log samtools          4
expect "samtools fastq, single-end"           samtools_single.log  samtools          298003
expect "bowtie2 --un-conc, 2 x pairs"         bowtie2_paired.log   bowtie2-paired    3695082
expect "bowtie2 --un"                         bowtie2_unpaired.log bowtie2-unpaired  288005
expect "bowtie2 --un-conc, empty input"       bowtie2_no_reads.log bowtie2-paired    0
expect "bowtie2 --un, empty input"            bowtie2_no_reads.log bowtie2-unpaired  0
# The paired summary's mate lines ("aligned 0 times") must not be taken
# for unpaired input, nor a missing line for zero
expect_error "samtools mode on a bowtie2 log"     bowtie2_paired.log   samtools
expect_error "bowtie2-paired mode on an unpaired log" bowtie2_unpaired.log bowtie2-paired
expect_error "empty log"                          empty.log            samtools
expect_error "unknown mode"                       samtools_paired.log  bwa

echo ""
echo "${PASS} passed, ${FAIL} failed"
[ "${FAIL}" -eq 0 ]
