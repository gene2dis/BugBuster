#!/bin/bash
#
# Regression tests for bin/gtdbtk_split.py, the per-sample split of
# GTDB_TK_BATCH's pooled GTDB-Tk summaries (roadmap Phase 20, RV3; review
# finding R4-02): every pooled bin is named <sample>__<bin>, and a bin whose
# prefix is not a sample of the batch, or that has no prefix, must stop the
# task instead of being dropped (a sample id holding '__' used to lose its rows
# or have them written into another sample's report). The script needs only
# the Python standard library, so it runs on the host python3.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
SPLIT="${REPO_DIR}/bin/gtdbtk_split.py"

TMP_DIR=$(mktemp -d)
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

PASS=0
FAIL=0
ok() { echo "✓ $1"; PASS=$((PASS + 1)); }
ko() { echo "✗ $1"; FAIL=$((FAIL + 1)); }

echo "=== gtdbtk_split.py test suite ==="

HEADER=$'user_genome\tclassification\tnote'
summary() {
    # summary <file> <user_genome>...: a pooled summary with one row per name
    local file="$1"; shift
    printf '%s\n' "${HEADER}" > "${file}"
    local name
    for name in "$@"; do
        printf '%s\td__Bacteria;p__Bacillota\tN/A\n' "${name}" >> "${file}"
    done
}

# A good split: rows go to their sample, the prefix is removed, a sample
# without bins gets a header-only table, and an id holding '_' (S_1) is kept
mkdir -p "${TMP_DIR}/good" && cd "${TMP_DIR}/good"
summary pooled.tsv "S_1__S_1_binette_bin1" "S2__SemiBin_3" "S_1__S_1_binette_bin2"
if python3 "${SPLIT}" --samples S_1 S2 S3 --domain bac120 --summary pooled.tsv > split.log 2>&1; then
    ok "good split exits 0"
else
    ko "good split exits 0: $(cat split.log)"
fi
[ "$(cut -f1 S_1_gtdbtk_bac120.tsv | tr '\n' ' ')" = "user_genome S_1_binette_bin1 S_1_binette_bin2 " ] \
    && ok "sample S_1 keeps its two bins, prefix removed" || ko "sample S_1 rows: $(cut -f1 S_1_gtdbtk_bac120.tsv | tr '\n' ' ')"
[ "$(sed -n 2p S2_gtdbtk_bac120.tsv)" = $'SemiBin_3\td__Bacteria;p__Bacillota\tN/A' ] \
    && ok "sample S2 row unchanged apart from the prefix" || ko "sample S2 row: $(sed -n 2p S2_gtdbtk_bac120.tsv)"
[ "$(cat S3_gtdbtk_bac120.tsv)" = "${HEADER}" ] && ok "sample without bins: header-only table" || ko "sample S3 table"

# The domain without a summary file: header-only tables in GTDB-Tk's columns
python3 "${SPLIT}" --samples S_1 S2 --domain ar53 --summary gtdbtk.ar53.summary.tsv > /dev/null
[ "$(head -c 24 S2_gtdbtk_ar53.tsv)" = $'user_genome\tclassificati' ] && [ "$(wc -l < S2_gtdbtk_ar53.tsv)" -eq 1 ] \
    && ok "missing domain summary: header-only tables" || ko "missing domain summary"

# A prefix that is not a sample of the batch (before, the rows were dropped
# with no message: a sample id 'site__rep1' parsed as 'site'). An id holding
# '__' cannot be split back unambiguously at all, which is why INPUT_CHECK
# rejects such ids
mkdir -p "${TMP_DIR}/foreign" && cd "${TMP_DIR}/foreign"
summary pooled.tsv "S1__bin1" "S9__bin1"
if python3 "${SPLIT}" --samples S1 S2 --domain bac120 --summary pooled.tsv > split.log 2>&1; then
    ko "unknown prefix S9 fails"
else
    grep -q "bin 'S9__bin1' has the prefix 'S9', which is not a sample of this batch (S1, S2)" split.log \
        && ok "unknown prefix S9 fails and names the bin" || ko "unknown prefix message: $(cat split.log)"
fi

# A bin name without the '<sample>__' prefix
mkdir -p "${TMP_DIR}/noprefix" && cd "${TMP_DIR}/noprefix"
summary pooled.tsv "S1__bin1" "bin2"
if python3 "${SPLIT}" --samples S1 --domain ar53 --summary pooled.tsv > split.log 2>&1; then
    ko "bin without prefix fails"
else
    grep -q "bin 'bin2' has no '<sample>__' prefix" split.log \
        && ok "bin without prefix fails and names the bin" || ko "no-prefix message: $(cat split.log)"
fi

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
