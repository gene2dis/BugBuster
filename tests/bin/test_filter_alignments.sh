#!/bin/bash
#
# Direct tests for bin/filter_alignments.sh, the identity / aligned-length
# filter between bwa-mem3 and samtools in BWAMEM3_SAMTOOLS. Asserts:
#   - header lines pass through unchanged,
#   - secondary and supplementary records are dropped,
#   - unmapped records pass through unchanged,
#   - a long, near-identical alignment is kept unchanged,
#   - a read shorter than 50 bp counts when >= 80% of it aligns (a 40 bp
#     read with 40M or 32M8S is kept, 31M9S is not),
#   - a short (< 50 M bases) or divergent (< 95% identity) alignment becomes
#     an unmapped record: flag 0x4 set, 0x2 cleared, MAPQ 0, CIGAR *, TLEN 0,
#     optional fields dropped, RNAME/POS/mate fields/SEQ/QUAL kept,
#   - insertions and deletions count in the identity denominator,
#   - thresholds 0 0 keep every primary alignment,
#   - out-of-range thresholds fail loudly.
#
# Runs on the host (bash + awk only; the container has mawk, so no gawk
# bit functions are allowed in the script).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
FILTER="${REPO_DIR}/bin/filter_alignments.sh"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0

check() {
    local desc="$1" cond="$2"
    if eval "${cond}"; then
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${desc}"
        FAIL=$((FAIL + 1))
    fi
}

# Record <name> <flag> <mapq> <cigar> <tlen> [tags...]: a 100 bp read on c1
# (READ_LEN overrides the length)
rec() {
    local name=$1 flag=$2 mapq=$3 cigar=$4 tlen=$5
    shift 5
    local seq qual len=${READ_LEN:-100}
    seq=$(printf 'A%.0s' $(seq 1 "${len}"))
    qual=$(printf 'I%.0s' $(seq 1 "${len}"))
    local line="${name}\t${flag}\tc1\t101\t${mapq}\t${cigar}\t=\t301\t${tlen}\t${seq}\t${qual}"
    for t in "$@"; do line="${line}\t${t}"; done
    printf '%b\n' "${line}"
}

{
    printf '@SQ\tSN:c1\tLN:1000\n'
    printf '@PG\tID:bwa-mem3\tPN:bwa-mem3\n'
    rec keep      99  60 100M        300 NM:i:2 AS:i:94
    rec short     99  60 40M60S      300 NM:i:0 AS:i:40
    rec divergent 99  60 100M        300 NM:i:8 AS:i:60
    rec indel     99  60 50M5I45M    300 NM:i:5 AS:i:70
    rec del_ok    99  60 60M3D40M    300 NM:i:3 AS:i:85
    rec edge      99  60 20S80M      300 NM:i:1 AS:i:75
    rec second    355 0  100M        300 NM:i:0
    rec supp      2147 60 30M70H     0   NM:i:0
    rec unmapped  69  0  '*'         0
    rec no_nm     99  60 100M        300 AS:i:100
} > "${TMP_DIR}/in.sam"

# 40 bp reads (2x50 library, or fastp-trimmed): 80% of 40 = 32 aligned bases
{
    printf '@SQ\tSN:c1\tLN:1000\n'
    READ_LEN=40 rec short_full  99 60 40M    100 NM:i:0 AS:i:40
    READ_LEN=40 rec short_80    99 60 32M8S  100 NM:i:1 AS:i:30
    READ_LEN=40 rec short_local 99 60 31M9S  100 NM:i:0 AS:i:31
    READ_LEN=40 rec short_div   99 60 40M    100 NM:i:3 AS:i:25
} > "${TMP_DIR}/short.sam"

echo "=== filter_alignments.sh test suite ==="

"${FILTER}" 0.95 50 < "${TMP_DIR}/in.sam" > "${TMP_DIR}/out.sam"
out="${TMP_DIR}/out.sam"
field() { awk -F'\t' -v n="$1" -v f="$2" '$1 == n { print $f }' "${out}"; }

check "header lines pass through" \
    "[[ \$(grep -c '^@' '${out}') -eq 2 ]]"
check "secondary record dropped" \
    "! grep -q '^second' '${out}'"
check "supplementary record dropped" \
    "! grep -q '^supp' '${out}'"
check "one record per remaining read" \
    "[[ \$(grep -vc '^@' '${out}') -eq 8 ]]"
check "unmapped record unchanged" \
    "[[ \"\$(grep '^unmapped' '${out}')\" == \"\$(grep '^unmapped' '${TMP_DIR}/in.sam')\" ]]"
check "100M with 98% identity kept unchanged" \
    "[[ \"\$(grep '^keep' '${out}')\" == \"\$(grep '^keep' '${TMP_DIR}/in.sam')\" ]]"
check "80 aligned bases clipped at a contig end kept" \
    "[[ \$(field edge 6) == 20S80M ]]"
check "40 aligned bases -> unmapped (flag 99 + 4 - 2 = 101)" \
    "[[ \$(field short 2) == 101 ]]"
check "92% identity -> unmapped" \
    "[[ \$(field divergent 2) == 101 ]]"
check "insertion counts in identity: 5 / 100 columns = 95% kept" \
    "[[ \$(field indel 6) == 50M5I45M ]]"
check "deletion counts in identity: 3 / 103 columns kept" \
    "[[ \$(field del_ok 6) == 60M3D40M ]]"
check "unmapped conversion: MAPQ 0, CIGAR *, TLEN 0" \
    "[[ \$(field short 5) == 0 && \$(field short 6) == '*' && \$(field short 9) == 0 ]]"
check "unmapped conversion keeps RNAME/POS and mate fields" \
    "[[ \$(field short 3) == c1 && \$(field short 4) == 101 && \$(field short 7) == '=' && \$(field short 8) == 301 ]]"
check "unmapped conversion keeps SEQ and QUAL, drops optional fields" \
    "[[ \$(field short 10 | wc -c) -eq 101 && \$(awk -F'\t' '\$1 == \"short\" { print NF }' '${out}') -eq 11 ]]"
check "missing NM tag treated as identity 1" \
    "[[ \$(field no_nm 2) == 99 ]]"

"${FILTER}" 0.95 50 < "${TMP_DIR}/short.sam" > "${TMP_DIR}/short_out.sam"
sfield() { awk -F'\t' -v n="$1" -v f="$2" '$1 == n { print $f }' "${TMP_DIR}/short_out.sam"; }
check "40 bp read aligned end to end kept (80% rule, not 50 bases)" \
    "[[ \$(sfield short_full 2) == 99 ]]"
check "40 bp read with 32 of 40 bases aligned kept" \
    "[[ \$(sfield short_80 2) == 99 ]]"
check "40 bp read with 31 of 40 bases aligned -> unmapped" \
    "[[ \$(sfield short_local 2) == 101 ]]"
check "40 bp read at 92.5% identity -> unmapped" \
    "[[ \$(sfield short_div 2) == 101 ]]"

"${FILTER}" 0 0 < "${TMP_DIR}/in.sam" > "${TMP_DIR}/off.sam"
check "thresholds 0 0 keep every primary alignment unchanged" \
    "[[ \$(awk -F'\t' '!/^@/ && int(\$2 / 4) % 2 == 0' '${TMP_DIR}/off.sam' | wc -l) -eq 7 ]]"

check "min_identity above 1 fails" \
    "! '${FILTER}' 1.5 50 < '${TMP_DIR}/in.sam' > /dev/null 2>&1"
check "negative min_aligned_length fails" \
    "! '${FILTER}' 0.95 -1 < '${TMP_DIR}/in.sam' > /dev/null 2>&1"

echo ""
echo "Passed: ${PASS}  Failed: ${FAIL}"
[[ ${FAIL} -eq 0 ]]
