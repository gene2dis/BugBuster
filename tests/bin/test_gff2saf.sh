#!/bin/bash
#
# Direct tests for bin/gff2saf.sh (functional annotation T3, design doc
# Section 4.4): the Prodigal/Pyrodigal GFF -> featureCounts SAF converter used
# by FEATURECOUNTS_GENES. Asserts:
#   - the GeneID rewrite from Prodigal's <seqnum>_<genenum> ID attribute to
#     the FAA-matching <contig>_<genenum> form (the T4 join precondition),
#   - gzipped and plain GFF input both accepted,
#   - an empty (header-only) GFF yields a header-only SAF,
#   - Chr/Start/End/Strand columns carried through verbatim,
#   - a CDS row with a malformed ID attribute fails loudly.
#
# Runs on the host (bash + awk + gzip only), against the committed fixtures
# in tests/data/gff/ (regenerate with tests/bin/make_gff_fixtures.sh).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
GFF2SAF="${REPO_DIR}/bin/gff2saf.sh"
FIXTURES="${REPO_DIR}/tests/data/gff"

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

echo "=== gff2saf.sh test suite ==="

#
# 1. Gzipped fixture: header + GeneID rewrite + column carry-through
#
"${GFF2SAF}" "${FIXTURES}/test_contig1_genes.gff.gz" "${TMP_DIR}/genes.saf"

check "gzipped GFF accepted, SAF has header + 2 data rows" \
    '[ "$(wc -l < "${TMP_DIR}/genes.saf")" -eq 3 ] && [ "$(head -1 "${TMP_DIR}/genes.saf")" = "$(printf "GeneID\tChr\tStart\tEnd\tStrand")" ]'

check "GeneID rewritten from ID=1_1 to contig_1_1" \
    '[ "$(sed -n 2p "${TMP_DIR}/genes.saf" | cut -f1)" = "contig_1_1" ]'

check "GeneID rewritten from ID=1_2 to contig_1_2" \
    '[ "$(sed -n 3p "${TMP_DIR}/genes.saf" | cut -f1)" = "contig_1_2" ]'

check "Chr/Start/End/Strand carried through (row 1: contig_1 1 300 +)" \
    '[ "$(sed -n 2p "${TMP_DIR}/genes.saf" | cut -f2-5)" = "$(printf "contig_1\t1\t300\t+")" ]'

check "Chr/Start/End/Strand carried through (row 2: contig_1 601 900 -)" \
    '[ "$(sed -n 3p "${TMP_DIR}/genes.saf" | cut -f2-5)" = "$(printf "contig_1\t601\t900\t-")" ]'

#
# 2. Plain (uncompressed) input gives identical output
#
gzip -cd "${FIXTURES}/test_contig1_genes.gff.gz" > "${TMP_DIR}/genes.gff"
"${GFF2SAF}" "${TMP_DIR}/genes.gff" "${TMP_DIR}/genes_plain.saf"

check "plain GFF input produces identical SAF" \
    'cmp -s "${TMP_DIR}/genes.saf" "${TMP_DIR}/genes_plain.saf"'

#
# 3. Empty (header-only) GFF -> header-only SAF
#
"${GFF2SAF}" "${FIXTURES}/test_empty_genes.gff.gz" "${TMP_DIR}/empty.saf"

check "empty GFF yields header-only SAF" \
    '[ "$(wc -l < "${TMP_DIR}/empty.saf")" -eq 1 ] && [ "$(head -1 "${TMP_DIR}/empty.saf")" = "$(printf "GeneID\tChr\tStart\tEnd\tStrand")" ]'

#
# 4. Malformed ID attribute fails loudly
#
printf '##gff-version  3\ncontig_1\tpyrodigal\tCDS\t1\t300\t1.0\t+\t0\tpartial=00;start_type=ATG\n' \
    > "${TMP_DIR}/bad.gff"

check "CDS row without an ID attribute exits non-zero" \
    '! "${GFF2SAF}" "${TMP_DIR}/bad.gff" "${TMP_DIR}/bad.saf" 2>/dev/null'

#
# Summary
#
echo ""
echo "Passed: ${PASS}  Failed: ${FAIL}"
[ "${FAIL}" -eq 0 ]
