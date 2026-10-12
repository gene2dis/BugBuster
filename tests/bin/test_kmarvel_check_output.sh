#!/bin/bash
#
# Regression tests for bin/kmarvel_check_output.sh, which turns two K-MARVEL
# 1.0.0 exit-0 outcomes into what the KMARVEL module needs:
#   - zero ARGs: kmarvel leaves a header-less, empty TSV -> the script writes
#     the header (it must match the header of a real K-MARVEL table)
#   - "Protein BLAST failed" with ARGs found -> exit 1 (unfiltered results)
#   - no "Found N Antimicrobial Resistance Genes" line -> exit 1
# Fixtures in tests/data/kmarvel/ are real K-MARVEL stdout logs and tables
# from the pipeline image (sampleA: 2 ARGs; sampleB: 0 ARGs).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
CHECK="${REPO_DIR}/bin/kmarvel_check_output.sh"
FIXTURES="${REPO_DIR}/tests/data/kmarvel"

TMP_DIR=$(mktemp -d)
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

PASS=0
FAIL=0

check() {
    # check <desc> <command...>
    local desc="$1"
    shift
    if "$@"; then
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${desc}"
        FAIL=$((FAIL + 1))
    fi
}

fails() {
    ! "$@" > /dev/null 2>&1
}

echo "=== kmarvel_check_output.sh test suite ==="
cd "${TMP_DIR}"

# ARGs found, BLAST ran: the table is left as kmarvel wrote it
cp "${FIXTURES}/sampleA_kmarvel.log" a.log
cp "${FIXTURES}/sampleA_kmarvel.tsv" a.tsv
check "ARGs found with BLAST columns: passes" "${CHECK}" a.log a.tsv
check "ARGs found: table left unchanged" cmp -s a.tsv "${FIXTURES}/sampleA_kmarvel.tsv"

# Zero ARGs: kmarvel's empty TSV gets the header of a real table
cp "${FIXTURES}/sampleB_kmarvel.log" b.log
: > b.tsv
check "zero ARGs: passes" "${CHECK}" b.log b.tsv
check "zero ARGs: header-only table written" test "$(wc -l < b.tsv)" -eq 1
check "zero ARGs: header matches a real K-MARVEL table" \
    test "$(cat b.tsv)" = "$(head -1 "${FIXTURES}/sampleA_kmarvel.tsv")"

# BLAST failed although ARGs were found: kmarvel exits 0, the task must not
sed 's/^\[INFO\] Processed .*/Protein BLAST failed with exit code: 2/' a.log > blastfail.log
cp "${FIXTURES}/sampleA_kmarvel.tsv" blastfail.tsv
check "BLAST failure with ARGs found: exits non-zero" fails "${CHECK}" blastfail.log blastfail.tsv
sed 's/^\[INFO\] Processed .*/Protein BLAST terminated by signal./' a.log > blastsig.log
check "BLAST killed with ARGs found: exits non-zero" fails "${CHECK}" blastsig.log blastfail.tsv

# ARGs found but the table lacks the BLAST columns (unfiltered layout)
cut -f1-14 "${FIXTURES}/sampleA_kmarvel.tsv" > noblast.tsv
check "ARGs found, no blast_e_value column: exits non-zero" fails "${CHECK}" a.log noblast.tsv

# kmarvel stopped before reporting its ARG count
grep -v 'Found .* Antimicrobial Resistance Genes' a.log > truncated.log
check "no ARG-count line: exits non-zero" fails "${CHECK}" truncated.log a.tsv

check "wrong argument count: exits non-zero" fails "${CHECK}" a.log

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
