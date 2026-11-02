#!/bin/bash
#
# Keeps docs/WORKFLOW_DIAGRAM.md in step with the code (roadmap D1): every process
# or subworkflow that main.nf or a local subworkflow includes (under its alias,
# e.g. NFCORE_SAMTOOLS_INDEX) must be named in one of the file's Mermaid
# diagrams or tables (prose does not count), and every upper-case name the file
# writes in backticks must still be included somewhere. Functions (flagOn,
# validateParameters, ...) are camelCase and are not checked.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
DOC="${DIAGRAM_DOC:-${REPO_DIR}/docs/WORKFLOW_DIAGRAM.md}"

TMP_DIR=$(mktemp -d)
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

PASS=0
FAIL=0

# Names included by main.nf and the local subworkflows: the alias after "as"
# when there is one, else the name; include blocks may span lines and hold
# several names separated by ";"
included_names() {
    cat "${REPO_DIR}/main.nf" "${REPO_DIR}"/subworkflows/local/*.nf \
        | tr '\n' ' ' \
        | grep -oE 'include[[:space:]]*\{[^}]*\}' \
        | sed -E 's/^include[[:space:]]*\{//; s/\}$//' \
        | tr ';' '\n' \
        | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//' \
        | awk 'NF { print ($2 == "as" ? $3 : $1) }' \
        | grep -E '^[A-Z][A-Z0-9_]*$' \
        | sort -u
}

# check_doc <doc>: prints one line per missing or stale name, returns 1 if any
check_doc() {
    local doc="$1" bad=0 name
    # the Mermaid blocks and the table rows
    awk '/^```mermaid/ { m = 1; next } m && /^```/ { m = 0 } m || /^\|/' "${doc}" > "${TMP_DIR}/named.txt"
    while read -r name; do
        if ! grep -qw -- "${name}" "${TMP_DIR}/named.txt"; then
            echo "missing: ${name} is included by the code but not in a diagram or table of ${doc##*/}"
            bad=1
        fi
    done < "${TMP_DIR}/included.txt"
    # backticked upper-case names are process / subworkflow names
    while read -r name; do
        if ! grep -qx -- "${name}" "${TMP_DIR}/included.txt"; then
            echo "stale: ${doc##*/} names \`${name}\`, which the code no longer includes"
            bad=1
        fi
    done < <(grep -oE '`[A-Z][A-Z0-9_]+`' "${doc}" | tr -d '`' | sort -u)
    return "${bad}"
}

included_names > "${TMP_DIR}/included.txt"
n_names=$(wc -l < "${TMP_DIR}/included.txt")

# The parser must find the stages main.nf is built from
if [ "${n_names}" -ge 80 ] && grep -qx RESULTS_SUMMARY "${TMP_DIR}/included.txt" \
        && grep -qx NFCORE_SAMTOOLS_INDEX "${TMP_DIR}/included.txt" \
        && grep -qx BWAMEM3_BUILD_COMBINED "${TMP_DIR}/included.txt" \
        && ! grep -qx SAMTOOLS_INDEX "${TMP_DIR}/included.txt"; then
    echo "✓ ${n_names} process and subworkflow names read from main.nf and subworkflows/local (aliases resolved)"
    PASS=$((PASS + 1))
else
    echo "✗ include parsing: ${n_names} names, expected >= 80 with RESULTS_SUMMARY, NFCORE_SAMTOOLS_INDEX, BWAMEM3_BUILD_COMBINED and no unaliased SAMTOOLS_INDEX"
    FAIL=$((FAIL + 1))
fi

# The real document
if check_doc "${DOC}" > "${TMP_DIR}/doc.log"; then
    echo "✓ ${DOC#"${REPO_DIR}"/} names every included process and subworkflow, and no stale one"
    PASS=$((PASS + 1))
else
    echo "✗ ${DOC#"${REPO_DIR}"/} is out of step with the code:"
    sed 's/^/    /' "${TMP_DIR}/doc.log"
    FAIL=$((FAIL + 1))
fi

# The check itself fails on a missing and on a stale name
sed -E 's/\bMICROBECENSUS\b/MICROBE_CENSUS_GONE/g' "${DOC}" > "${TMP_DIR}/missing.md"
if ! check_doc "${TMP_DIR}/missing.md" > "${TMP_DIR}/missing.log" && grep -q 'missing: MICROBECENSUS' "${TMP_DIR}/missing.log"; then
    echo "✓ a process left out of the document is reported"
    PASS=$((PASS + 1))
else
    echo "✗ removing MICROBECENSUS from the document was not reported"
    FAIL=$((FAIL + 1))
fi
for stale in PRODIGAL_CONTIGS KARGA; do
    { cat "${DOC}"; echo "Old step: \`${stale}\`"; } > "${TMP_DIR}/stale.md"
    if ! check_doc "${TMP_DIR}/stale.md" > "${TMP_DIR}/stale.log" && grep -q "stale: .*${stale}" "${TMP_DIR}/stale.log"; then
        echo "✓ ${stale}, which the code no longer includes, is reported"
        PASS=$((PASS + 1))
    else
        echo "✗ a stale \`${stale}\` was not reported"
        FAIL=$((FAIL + 1))
    fi
done
# A name left only in prose (QFILTER taken out of every diagram and table, then
# mentioned in a sentence) counts as missing
sed -E '/^\|/ s/\bQFILTER\b/QF/g; /^```mermaid/,/^```/ s/\bQFILTER\b/QF/g' "${DOC}" > "${TMP_DIR}/prose.md"
echo 'QFILTER is mentioned in prose only.' >> "${TMP_DIR}/prose.md"
if ! check_doc "${TMP_DIR}/prose.md" > "${TMP_DIR}/prose.log" && grep -q 'missing: QFILTER' "${TMP_DIR}/prose.log"; then
    echo "✓ a process named only in prose is reported"
    PASS=$((PASS + 1))
else
    echo "✗ QFILTER named only in prose was not reported"
    FAIL=$((FAIL + 1))
fi

echo
echo "${PASS} passed, ${FAIL} failed"
[ "${FAIL}" -eq 0 ]
