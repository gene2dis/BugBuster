#!/bin/bash
#
# Regression test for design doc Q21 (T9 finding, 2026-10-09): under
# singularity/apptainer the BLOBTOOLS image is read-only, and `blobtools
# create --nodes --names` wrote its parsed nodesDB into the blobtools package
# directory by default (OSError: Read-only file system). The module now passes
# `--db nodesDB.txt` so the nodesDB is written in the task directory.
#
# The pinned image runs under `docker run --read-only` (a read-only root, as in
# a sif) on a tiny on-the-fly fixture:
#   - BLOBTOOLS_NODESDB builds nodesDB.txt in the task directory once, with
#     BlobTools' own parser (roadmap S1), and BLOBTOOLS reads it with --db;
#     the parse runs on a read-only image and create + view work from it,
#   - create (with --db) + view succeed and write nodesDB.txt in the work dir,
#   - create without --db fails with the read-only error (the test can see
#     the defect),
#   - with a 2025-style taxdump ('domain' top rank, Viruses an 'acellular
#     root' with a realm below it) the superkingdom column is 'undef' as is,
#     and Bacteria / Viruses after bin/blobtools_nodes_fix.sh (which the
#     BLOBTOOLS_NODESDB runs).
#
# Requirements: docker (present on GitHub ubuntu-latest runners).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
MODULE="${REPO_DIR}/modules/local/blobtools/main.nf"
NODESDB_MODULE="${REPO_DIR}/modules/local/blobtools_nodesdb/main.nf"

# Use the exact image pinned in the module so the test cannot drift from it
# (the docker URI of the engine closure; the depot sif is the same build)
BLOB_IMG=$(grep -o "'quay.io/biocontainers/blobtools:[^']*'" "${MODULE}" | tr -d "'" | head -1)
if [ -z "${BLOB_IMG}" ]; then
    echo "ERROR: could not extract the pinned BLOBTOOLS image from ${MODULE}"
    exit 1
fi

echo "=== BLOBTOOLS read-only image test suite ==="
echo "Image: ${BLOB_IMG}"
echo ""

TMP_DIR=$(mktemp -d)
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

PASS=0
FAIL=0
pass() { echo "PASS: $1"; PASS=$((PASS + 1)); }
fail() { echo "FAIL: $1"; FAIL=$((FAIL + 1)); }

run_ro() {
    # run_ro <workdir> <log> <blobtools args...> — read-only root, writable /tmp
    local workdir="$1" log="$2"
    shift 2
    docker run --rm --read-only --tmpfs /tmp -u "$(id -u):$(id -g)" -e HOME=/tmp \
        -v "${workdir}:/w" -w /w "${BLOB_IMG}" blobtools "$@" > "${log}" 2>&1
}

make_fixture() {
    local d="$1"
    mkdir -p "${d}"
    printf ">c1\nACGTACGTACGTACGTACGTAAACCCGGGTTT\n>c2\nGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT\n" > "${d}/contigs.fa"
    printf "c1\t562\t200.0\nc2\t1280\t150.0\n" > "${d}/hits.tsv"
    printf "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n562\t|\t2\t|\tspecies\t|\n1280\t|\t2\t|\tspecies\t|\n" > "${d}/nodes.dmp"
    printf "1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\t\t|\tscientific name\t|\n562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n1280\t|\tStaphylococcus aureus\t|\t\t|\tscientific name\t|\n" > "${d}/names.dmp"
    printf "## blobtools v1.1.1\n## Total Reads = 20\n# contig_id\tread_cov\tbase_cov\nc1\t10\t9.5\nc2\t10\t9.0\n" > "${d}/reads.cov"
}

#
# 1. The module keeps the fix
#
if grep -q -- '--db ${nodesdb}' "${MODULE}" && ! grep -q -- '--nodes' "${MODULE}"; then
    pass "module: BLOBTOOLS reads the prebuilt nodesDB with --db, no --nodes"
else
    fail "module: BLOBTOOLS reads the prebuilt nodesDB with --db, no --nodes"
fi
if grep -q "nodesDB='nodesDB.txt'" "${NODESDB_MODULE}" \
        && grep -q "'quay.io/biocontainers/blobtools:1.1.1--py_1'" "${NODESDB_MODULE}"; then
    pass "module: BLOBTOOLS_NODESDB writes nodesDB.txt in the task dir, same image"
else
    fail "module: BLOBTOOLS_NODESDB writes nodesDB.txt in the task dir, same image"
fi

#
# 1b. The module's nodesDB parse on a read-only image, then create --db from it
#
P="${TMP_DIR}/prebuilt"
make_fixture "${P}"
cp "${P}/nodes.dmp" "${P}/nodes_blobtools.dmp"
# the python line of BLOBTOOLS_NODESDB, with the task's names.dmp path
PARSE=$(grep -o 'import blobtools.BtIO.*nodesDBdefault=None)' "${NODESDB_MODULE}" | sed 's/\\\${names_dmp}/names.dmp/')
if docker run --rm --read-only --tmpfs /tmp -u "$(id -u):$(id -g)" -e HOME=/tmp \
        -v "${P}:/w" -w /w "${BLOB_IMG}" python -c "${PARSE}" > "${P}/parse.log" 2>&1 \
        && [ -s "${P}/nodesDB.txt" ]; then
    pass "BLOBTOOLS_NODESDB parse writes nodesDB.txt on a read-only image"
else
    fail "BLOBTOOLS_NODESDB parse writes nodesDB.txt on a read-only image"; cat "${P}/parse.log"
fi
if run_ro "${P}" "${P}/create.log" create --infile contigs.fa --hitsfile hits.tsv \
        --db nodesDB.txt --cov reads.cov --out p \
        && run_ro "${P}" "${P}/view.log" view --input p.blobDB.json --out p_Blob_table --rank all \
        && grep -q "Escherichia coli" "${P}"/p_Blob_table*; then
    pass "create --db with the prebuilt nodesDB alone gives the species call"
else
    fail "create --db with the prebuilt nodesDB alone gives the species call"; cat "${P}/create.log" "${P}/view.log" 2>/dev/null
fi

#
# 2. With --db: create + view succeed on a read-only root
#
W="${TMP_DIR}/with_db"
make_fixture "${W}"
if run_ro "${W}" "${W}/create.log" create --infile contigs.fa --hitsfile hits.tsv \
        --nodes nodes.dmp --names names.dmp --db nodesDB.txt --cov reads.cov --out t; then
    pass "create with --db succeeds on a read-only image"
else
    fail "create with --db succeeds on a read-only image"; cat "${W}/create.log"
fi
if [ -s "${W}/nodesDB.txt" ]; then
    pass "create with --db writes nodesDB.txt in the work dir"
else
    fail "create with --db writes nodesDB.txt in the work dir"
fi
if run_ro "${W}" "${W}/view.log" view --input t.blobDB.json --out t_Blob_table --rank all \
        && grep -q "Escherichia coli" "${W}"/t_Blob_table*; then
    pass "view writes a blob table with the species call"
else
    fail "view writes a blob table with the species call"; cat "${W}/view.log"
fi

#
# 3. Without --db: the read-only failure is visible (negative control)
#
N="${TMP_DIR}/without_db"
make_fixture "${N}"
if run_ro "${N}" "${N}/create.log" create --infile contigs.fa --hitsfile hits.tsv \
        --nodes nodes.dmp --names names.dmp --cov reads.cov --out u; then
    fail "create without --db fails on a read-only image (negative control)"
elif grep -q "Read-only file system" "${N}/create.log"; then
    pass "create without --db fails on a read-only image (negative control)"
else
    fail "create without --db fails with the read-only error"; cat "${N}/create.log"
fi

#
# 4. 2025 NCBI taxdump: 'domain' and the Viruses 'acellular root'
#
make_fixture_2025() {
    local d="$1"
    make_fixture "${d}"
    printf "c1\t562\t200.0\nc2\t10760\t150.0\n" > "${d}/hits.tsv"
    printf "1\t|\t1\t|\tno rank\t|\n131567\t|\t1\t|\tcellular root\t|\n2\t|\t131567\t|\tdomain\t|\n562\t|\t2\t|\tspecies\t|\n10239\t|\t1\t|\tacellular root\t|\n2731341\t|\t10239\t|\trealm\t|\n10760\t|\t2731341\t|\tspecies\t|\n" > "${d}/nodes.dmp"
    printf "1\t|\troot\t|\t\t|\tscientific name\t|\n131567\t|\tcellular organisms\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\t\t|\tscientific name\t|\n562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n10239\t|\tViruses\t|\t\t|\tscientific name\t|\n2731341\t|\tDuplodnaviria\t|\t\t|\tscientific name\t|\n10760\t|\tEscherichia phage T7\t|\t\t|\tscientific name\t|\n" > "${d}/names.dmp"
}
superkingdoms() {
    # superkingdoms <blob table> — c1,c2 superkingdom.t column values
    grep -v '^#' "$1" | sort | cut -f6 | paste -sd, -
}
if grep -q "blobtools_nodes_fix.sh" "${NODESDB_MODULE}"; then
    pass "module: nodes.dmp passes through blobtools_nodes_fix.sh"
else
    fail "module: nodes.dmp passes through blobtools_nodes_fix.sh"
fi
R="${TMP_DIR}/taxdump_2025_raw"
make_fixture_2025 "${R}"
run_ro "${R}" "${R}/create.log" create --infile contigs.fa --hitsfile hits.tsv \
    --nodes nodes.dmp --names names.dmp --db nodesDB.txt --cov reads.cov --out r \
    && run_ro "${R}" "${R}/view.log" view --input r.blobDB.json --out r_Blob_table --rank all || true
if [ "$(superkingdoms "${R}"/r_Blob_table*)" = "undef,undef" ]; then
    pass "2025 taxdump as is: superkingdom undef (the defect is visible)"
else
    fail "2025 taxdump as is: superkingdom undef (got $(superkingdoms "${R}"/r_Blob_table* 2>/dev/null))"
fi
F="${TMP_DIR}/taxdump_2025_fixed"
make_fixture_2025 "${F}"
"${REPO_DIR}/bin/blobtools_nodes_fix.sh" "${F}/nodes.dmp" > "${F}/nodes_blobtools.dmp"
if run_ro "${F}" "${F}/create.log" create --infile contigs.fa --hitsfile hits.tsv \
        --nodes nodes_blobtools.dmp --names names.dmp --db nodesDB.txt --cov reads.cov --out f \
        && run_ro "${F}" "${F}/view.log" view --input f.blobDB.json --out f_Blob_table --rank all \
        && [ "$(superkingdoms "${F}"/f_Blob_table*)" = "Bacteria,Viruses" ]; then
    pass "2025 taxdump after blobtools_nodes_fix.sh: superkingdom Bacteria, Viruses"
else
    fail "2025 taxdump after blobtools_nodes_fix.sh: superkingdom Bacteria, Viruses (got $(superkingdoms "${F}"/f_Blob_table* 2>/dev/null))"
fi
if cmp -s <("${REPO_DIR}/bin/blobtools_nodes_fix.sh" "${W}/nodes.dmp") "${W}/nodes.dmp"; then
    pass "blobtools_nodes_fix.sh leaves a 'superkingdom' taxdump unchanged"
else
    fail "blobtools_nodes_fix.sh leaves a 'superkingdom' taxdump unchanged"
fi

echo ""
echo "=== ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
