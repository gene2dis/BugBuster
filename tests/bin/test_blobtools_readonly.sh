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
#   - the module's BLOBTOOLS script still carries `--db nodesDB.txt`,
#   - create (with --db) + view succeed and write nodesDB.txt in the work dir,
#   - create without --db fails with the read-only error (the test can see
#     the defect).
#
# Requirements: docker (present on GitHub ubuntu-latest runners).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
MODULE="${REPO_DIR}/modules/local/blobtools/main.nf"

# Use the exact image pinned in the module so the test cannot drift from it
BLOB_IMG=$(grep -o "container '[^']*'" "${MODULE}" | head -1 | sed "s/container '//; s/'$//")
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
if grep -q -- "--db nodesDB.txt" "${MODULE}"; then
    pass "module: blobtools create passes --db nodesDB.txt"
else
    fail "module: blobtools create passes --db nodesDB.txt"
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

echo ""
echo "=== ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
