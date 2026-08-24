#!/bin/bash
#
# Regression tests for audit item #8: RGI_BWT's low-coverage fallback
# (bin/rgi_bwt_lowcov_fallback.sh) must create every output the module
# declares as non-optional — previously it skipped *.reference_mapping_stats.txt
# and *.sorted.length_100.bam, so the fallback itself caused a
# MissingFileException in the exact case it was written for. The bam emit is
# now optional (no alignments -> RGI_KMER skipped), and this test derives the
# expected file set from the module's own output block so a future
# non-optional output added without a fallback counterpart fails here.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
MODULE="${REPO_DIR}/modules/local/rgi_bwt/main.nf"
FALLBACK="${REPO_DIR}/bin/rgi_bwt_lowcov_fallback.sh"

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

echo "=== rgi_bwt_lowcov_fallback.sh test suite ==="

PREFIX="sampleX_rgi_bwt"
cd "${TMP_DIR}"
"${FALLBACK}" "${PREFIX}"

# Every non-optional path("*...") glob in the module's output block must be
# satisfied by a file the fallback created (versions.yml is written by the
# module itself; optional emits — e.g. the bam — are excluded).
mapfile -t globs < <(
    sed -n '/^    output:/,/^    when:/p' "${MODULE}" \
    | grep -v 'optional: *true' \
    | grep -o 'path("\*[^"]*")' \
    | sed 's/^path("\(.*\)")$/\1/'
)

check "module output block parsed (found some non-optional globs)" \
    test "${#globs[@]}" -ge 4

for glob in "${globs[@]}"; do
    matches=$(compgen -G "${glob}" || true)
    check "fallback creates a file matching non-optional output '${glob}'" \
        test -n "${matches}"
done

# The previously-missing file, asserted explicitly
check "reference_mapping_stats.txt exists and is non-empty" \
    test -s "${PREFIX}.reference_mapping_stats.txt"

# No BAM: the bam emit is optional and RGI_KMER must be skipped, never fed a fake
check "no fake BAM is fabricated" \
    bash -c "! compgen -G '*.bam' > /dev/null"

# The NA tables must stay parseable by RGI_REPORT (real tabs, NA row width == header width)
for table in allele_mapping_data gene_mapping_data; do
    f="${PREFIX}.${table}.txt"
    header_cols=$(head -1 "${f}" | awk -F'\t' '{print NF}')
    na_cols=$(sed -n '2p' "${f}" | awk -F'\t' '{print NF}')
    check "${table}: NA row width (${na_cols}) matches header width (${header_cols})" \
        test "${header_cols}" -eq "${na_cols}" -a "${header_cols}" -gt 1
done
check "allele table has the 27-column RGI 6.x header" \
    test "$(head -1 "${PREFIX}.allele_mapping_data.txt" | awk -F'\t' '{print NF}')" -eq 27
check "gene table has the 16-column RGI 6.x header" \
    test "$(head -1 "${PREFIX}.gene_mapping_data.txt" | awk -F'\t' '{print NF}')" -eq 16

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
