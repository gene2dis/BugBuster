#!/bin/bash
#
# Regression tests for audit item #6: CHECKM2_BATCH/GTDB_TK_BATCH must recognize
# every binner's output directory, with the MetaWRAP directory suffix derived
# from params.metawrap_completeness/contamination instead of a hardcoded
# "_metawrap_50_10_bins" literal. The recognition logic lives in
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
    # expect_id <desc> <expected "sample<TAB>binner"> <dir> <metawrap_suffix>
    local desc="$1" expected="$2" dir="$3" suffix="$4"
    local got
    if got=$("${IDENTIFY}" "${dir}" "${suffix}" 2>"${TMP_DIR}/stderr.log"); then
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
    # expect_reject <desc> <dir> <metawrap_suffix>
    local desc="$1" dir="$2" suffix="$3"
    if "${IDENTIFY}" "${dir}" "${suffix}" >/dev/null 2>&1; then
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
expect_id "metabat dir" "S1${TAB}metabat" "bin_1/S1_metabat_bins" "metawrap_50_10_bins"
expect_id "semibin dir (default --binners semibin, the audit-#6 case)" \
    "S1${TAB}semibin" "bin_1/S1_semibin_output_bins" "metawrap_50_10_bins"
expect_id "comebin dir" "S1${TAB}comebin" "bin_1/S1_comebin_bins" "metawrap_50_10_bins"

# Sample ids containing underscores and binner-like words
expect_id "underscore-heavy sample id" \
    "S10_Ago2021${TAB}semibin" "S10_Ago2021_semibin_output_bins" "metawrap_50_10_bins"
expect_id "sample id containing a binner keyword" \
    "metabat_test${TAB}metabat" "metabat_test_metabat_bins" "metawrap_50_10_bins"

# MetaWRAP with the default and non-default params.metawrap_* suffixes
expect_id "metawrap dir, default 50_10 suffix" \
    "S1${TAB}metawrap" "S1_metawrap_50_10_bins" "metawrap_50_10_bins"
expect_id "metawrap dir, custom 70_5 suffix (previously hardcoded to 50_10)" \
    "S10_Ago2021${TAB}metawrap" "S10_Ago2021_metawrap_70_5_bins" "metawrap_70_5_bins"

# COMEBin nested dir staged by its basename: sample id comes from the resolved path
mkdir -p "${TMP_DIR}/work/S2_comebin_bins/comebin_res/comebin_res_bins"
mkdir -p "${TMP_DIR}/staged"
ln -sf "${TMP_DIR}/work/S2_comebin_bins/comebin_res/comebin_res_bins" "${TMP_DIR}/staged/comebin_res_bins"
expect_id "comebin_res_bins symlink resolves sample from target path" \
    "S2${TAB}comebin" "${TMP_DIR}/staged/comebin_res_bins" "metawrap_50_10_bins"

# Failure modes: naming drift must be loud, not a silent skip
expect_reject "unrecognized directory fails" "S1_random_bins" "metawrap_50_10_bins"
expect_reject "metawrap dir with mismatched suffix fails" \
    "S1_metawrap_70_5_bins" "metawrap_50_10_bins"
mkdir -p "${TMP_DIR}/staged_bad"
ln -sf "${TMP_DIR}/work" "${TMP_DIR}/staged_bad/comebin_res_bins"
expect_reject "comebin_res_bins without a _comebin path component fails" \
    "${TMP_DIR}/staged_bad/comebin_res_bins" "metawrap_50_10_bins"

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
