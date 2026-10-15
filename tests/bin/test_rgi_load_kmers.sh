#!/bin/bash
#
# Regression tests for bin/rgi_load_kmers.sh (RGI_LOAD, RGI_LOAD_WILDCARD).
# rgi kmer_query (RGI_KMER, pathogen of origin) reads the k-mer database from
# localDB/, but neither database step loaded it: kmer_query logged an error,
# exited 0 and RGI_KMER wrote an empty '{}' for every sample. Runs the real
# `rgi load` inside the RGI image the modules pin, on tiny stand-in files
# (rgi load copies them without parsing), and checks:
#   - plain and gzipped k-mer files are loaded under kmer_query's names;
#   - the CARD/WildCARD entries of an earlier load are kept;
#   - a k-mer size the WildCARD directory lacks fails with a clear error.
#
# Requirements: docker.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
RGI_IMG=$(grep -o "'quay.io/biocontainers/rgi:[^']*'" "${REPO_DIR}/modules/local/rgi_load_wildcard/main.nf" \
    | tr -d "'" | head -1)

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

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
run_rgi() {
    # run_rgi <command...> — in ${RGI_IMG}, with the repo bin/ on PATH, in ${TMP_DIR}
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" -v "${TMP_DIR}:/work" -w /work \
        "${RGI_IMG}" bash -c 'export PATH="$PATH:/pipeline_bin"; "$@"' _ "$@"
}

echo "=== rgi_load_kmers.sh test suite (${RGI_IMG}) ==="

cd "${TMP_DIR}"
mkdir -p wildcard localDB
echo '{"_version": "9.9.9", "_comment": "test", "_timestamp": "2026-01-01"}' > card.json
echo '{"k": 1}' > wildcard/61_kmer_db.json
printf 'ACGT\t1\n' | gzip > wildcard/all_amr_61mers.txt.gz
# an earlier CARD + WildCARD load, as RGI_LOAD / RGI_LOAD_WILDCARD leave it
cat > localDB/loaded_databases.json <<'JSON'
{"card_json": {"data_version": "9.9.9", "model_type_used": ["homolog"]}, "card_canonical": {"data_version": "9.9.9", "model_type_used": ["homolog"]}, "card_variants": {"data_version": "9.9.9", "model_type_used": ["homolog"]}, "card_kmers": {"kmer_sizes": []}}
JSON

if run_rgi rgi_load_kmers.sh card.json wildcard 61 > load61.log 2>&1; then
    check "61-mer load succeeds" true
else
    check "61-mer load succeeds" false
    tail -5 load61.log
fi
check "k-mer database loaded as localDB/61mer_database.json" test -s localDB/61mer_database.json
check "gzipped AMR k-mers loaded as localDB/amr_61mer.txt" \
    bash -c "test \"\$(cat localDB/amr_61mer.txt)\" = \"\$(printf 'ACGT\t1')\""
check "loaded_databases.json lists k=61" grep -q '"kmer_sizes": \["61"\]' localDB/loaded_databases.json
check "earlier WildCARD load kept" grep -q '"card_variants": {"data_version": "9.9.9"' localDB/loaded_databases.json

if run_rgi rgi_load_kmers.sh card.json wildcard 31 > load31.log 2>&1; then
    check "missing 31-mer database fails" false
else
    check "missing 31-mer database fails" true
fi
check "failure names the missing file" grep -q "wildcard/31_kmer_db.json not found" load31.log

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
