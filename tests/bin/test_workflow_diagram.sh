#!/bin/bash
#
# Keeps docs/WORKFLOW_DIAGRAM.md in step with the code (roadmap D1): every process
# or subworkflow that main.nf or a local subworkflow includes (under its alias,
# e.g. NFCORE_SAMTOOLS_INDEX) must be named in the file, and every process name
# the file writes in backticks must still be included somewhere. Functions
# (flagOn, validateParameters, ...) are camelCase and are not checked.
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
    while read -r name; do
        if ! grep -qw -- "${name}" "${doc}"; then
            echo "missing: ${name} is included by the code but not named in ${doc##*/}"
            bad=1
        fi
    done < "${TMP_DIR}/included.txt"
    # backticked upper-case names with an underscore are process / subworkflow names
    while read -r name; do
        if ! grep -qx -- "${name}" "${TMP_DIR}/included.txt"; then
            echo "stale: ${doc##*/} names \`${name}\`, which the code no longer includes"
            bad=1
        fi
    done < <(grep -oE '`[A-Z][A-Z0-9]*_[A-Z0-9_]*`' "${doc}" | tr -d '`' | sort -u)
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
{ cat "${DOC}"; echo 'Old step: `PRODIGAL_CONTIGS`'; } > "${TMP_DIR}/stale.md"
if ! check_doc "${TMP_DIR}/stale.md" > "${TMP_DIR}/stale.log" && grep -q 'stale: .*PRODIGAL_CONTIGS' "${TMP_DIR}/stale.log"; then
    echo "✓ a process the code no longer includes is reported"
    PASS=$((PASS + 1))
else
    echo "✗ a stale \`PRODIGAL_CONTIGS\` was not reported"
    FAIL=$((FAIL + 1))
fi

echo
echo "${PASS} passed, ${FAIL} failed"
[ "${FAIL}" -eq 0 ]
