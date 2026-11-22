#!/bin/bash
#
# Regression tests for audit item #6: CHECKM2_BATCH/GTDB_TK_BATCH must recognize
# every binner's output directory and Binette's refined bins. The recognition
# logic lives in
# bin/identify_bin_dir.sh (on PATH inside the processes) so it can be tested
# directly, without containers or databases.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
IDENTIFY="${REPO_DIR}/bin/identify_bin_dir.sh"

TMP_DIR=$(mktemp -d)
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

PASS=0
FAIL=0

expect_id() {
    # expect_id <desc> <expected "sample<TAB>binner"> <dir> [sample_id]
    local desc="$1" expected="$2" dir="$3"
    shift 3
    local got
    if got=$("${IDENTIFY}" "${dir}" "$@" 2>"${TMP_DIR}/stderr.log"); then
        if [ "${got}" = "${expected}" ]; then
            echo "✓ ${desc}"
            PASS=$((PASS + 1))
        else
            echo "✗ ${desc} — expected '${expected}', got '${got}'"
            FAIL=$((FAIL + 1))
        fi
    else
        echo "✗ ${desc} — expected success, script failed:"
        sed 's/^/    /' "${TMP_DIR}/stderr.log"
        FAIL=$((FAIL + 1))
    fi
}

expect_reject() {
    # expect_reject <desc> <dir> [sample_id]
    local desc="$1" dir="$2"
    shift 2
    if "${IDENTIFY}" "${dir}" "$@" >/dev/null 2>&1; then
        echo "✗ ${desc} — expected nonzero exit, script succeeded"
        FAIL=$((FAIL + 1))
    else
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    fi
}

echo "=== identify_bin_dir.sh test suite ==="

TAB=$'\t'

# Plain binner directories (the shapes each binner module emits)
expect_id "metabat dir" "S1${TAB}metabat" "bin_1/S1_metabat_bins"
expect_id "semibin dir (default --binners semibin, the audit-#6 case)" \
    "S1${TAB}semibin" "bin_1/S1_semibin_output_bins"
expect_id "comebin dir" "S1${TAB}comebin" "bin_1/S1_comebin_bins"

# Sample ids containing underscores and binner-like words
expect_id "underscore-heavy sample id" \
    "site_A_rep1${TAB}semibin" "site_A_rep1_semibin_output_bins"
expect_id "sample id containing a binner keyword" \
    "metabat_test${TAB}metabat" "metabat_test_metabat_bins"

# Binette's refined bins
expect_id "binette dir" "S1${TAB}binette" "S1_binette_bins"
expect_id "binette dir, underscore-heavy sample id" \
    "site_A_rep1${TAB}binette" "site_A_rep1_binette_bins"

# The sample id passed by the process (roadmap Phase 20, RV5): COMEBin's nested
# comebin_res_bins dir takes it from the argument, however it was staged (a
# copy has no symlink to resolve; a work path holding '_comebin' gave its own
# prefix before); a named dir must carry the given id
mkdir -p "${TMP_DIR}/copied/comebin_res_bins"
expect_id "comebin_res_bins (copied, not a symlink) takes the given id" \
    "S2${TAB}comebin" "${TMP_DIR}/copied/comebin_res_bins" S2
mkdir -p "${TMP_DIR}/proj_comebin_test/work/S3_comebin_bins/comebin_res/comebin_res_bins"
expect_id "comebin_res_bins under a path holding _comebin takes the given id" \
    "S3${TAB}comebin" "${TMP_DIR}/proj_comebin_test/work/S3_comebin_bins/comebin_res/comebin_res_bins" S3
expect_id "named dir with its own id given" "S1${TAB}metabat" "bin_1/S1_metabat_bins" S1
expect_id "binette dir with its own id given" "site_A_rep1${TAB}binette" "site_A_rep1_binette_bins" site_A_rep1

# Failure modes: naming drift must be loud, not a silent skip
expect_reject "unrecognized directory fails" "S1_random_bins"
expect_reject "former MetaWRAP dir fails" "S1_metawrap_50_10_bins"
expect_reject "comebin_res_bins without a given id fails" "${TMP_DIR}/copied/comebin_res_bins"
expect_reject "named dir staged for another sample fails" "bin_1/S1_metabat_bins" S2

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
