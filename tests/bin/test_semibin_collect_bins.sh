#!/bin/bash
#
# Direct tests for bin/semibin_collect_bins.sh, the step after SemiBin2 that
# both SEMIBIN script branches call (roadmap Phase 29, CL3; issue #30): it
# moves SemiBin2's bins to the module output, writes SKIPPED.txt when the log
# shows SemiBin2 found no bins, and otherwise exits with SemiBin2's status.
# The log lines are SemiBin2 2.6.0's own messages (captured from the image). Asserts:
#   - status 0: the bin directory is moved to the output, bins kept;
#   - "Number of bins prior to reclustering: 0": SKIPPED.txt, exit 0, warning;
#   - "no must-link pairs can be generated" and "but all are shorter than 2500
#     basepairs" (SemiBin2's first step, no contig of 4,000 / 2,500 bp):
#     SKIPPED.txt, exit 0;
#   - any other failure: the same exit status, no output directory;
#   - a wrong argument count fails with a usage message.
#
# Runs on the host (bash + grep only); the SEMIBIN module nf-test runs the
# script inside the SemiBin2 image.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
COLLECT="${REPO_DIR}/bin/semibin_collect_bins.sh"

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

# run <case> <status> <log text>: in its own folder, with a SemiBin2-like
# output/output_bins holding one bin; sets RC and leaves stderr in err.txt
run() {
    local dir="${TMP_DIR}/$1"
    mkdir -p "${dir}/output/output_bins"
    echo ">c1" > "${dir}/output/output_bins/SemiBin_0.fa"
    printf '%s\n' "$3" > "${dir}/semibin.log"
    RC=0
    (cd "${dir}" && "${COLLECT}" "$2" output/output_bins s1_semibin_output_bins semibin.log s1) 2> "${dir}/err.txt" || RC=$?
    CASE="${dir}"
}

echo "=== semibin_collect_bins.sh test suite ==="

run success 0 "Binning finished"
check "status 0: exit 0" '[ "${RC}" -eq 0 ]'
check "status 0: the bins are moved to the output directory" '[ -f "${CASE}/s1_semibin_output_bins/SemiBin_0.fa" ] && [ ! -e "${CASE}/output/output_bins" ]'
check "status 0: no SKIPPED.txt" '[ ! -e "${CASE}/s1_semibin_output_bins/SKIPPED.txt" ]'

run zero_bins 1 "2026-11-27 12:00:00,000 - Number of bins prior to reclustering: 0"
check "zero bins: exit 0" '[ "${RC}" -eq 0 ]'
check "zero bins: SKIPPED.txt says SemiBin found no bins" '[ "$(cat "${CASE}/s1_semibin_output_bins/SKIPPED.txt")" = "SemiBin found no bins" ]'
check "zero bins: the output holds only SKIPPED.txt" '[ "$(ls "${CASE}/s1_semibin_output_bins")" = "SKIPPED.txt" ]'
check "zero bins: a warning names the sample" 'grep -q "WARNING: SemiBin found no bins for s1." "${CASE}/err.txt"'

run no_must_link 1 "2026-11-27 20:59:09 host SemiBin2[9] ERROR There are 4 contigs in input file contigs.fa, but none contain at least 4000 basepairs (the longest is 3400), so no must-link pairs can be generated."
check "no must-link pairs: exit 0 and SKIPPED.txt" '[ "${RC}" -eq 0 ] && [ -f "${CASE}/s1_semibin_output_bins/SKIPPED.txt" ]'

run all_short 1 "2026-11-27 20:58:45 host SemiBin2[10] ERROR Input file contigs.fa contains 3 contigs, but all are shorter than 2500 basepairs."
check "every contig under 2,500 bp: exit 0 and SKIPPED.txt" '[ "${RC}" -eq 0 ] && [ -f "${CASE}/s1_semibin_output_bins/SKIPPED.txt" ]'

run other_failure 137 "Killed"
check "other failure: SemiBin2's exit status (137) is kept" '[ "${RC}" -eq 137 ]'
check "other failure: no output directory" '[ ! -e "${CASE}/s1_semibin_output_bins" ]'

run other_failure_1 1 "RuntimeError: CUDA error"
check "other failure: exit status 1 is kept" '[ "${RC}" -eq 1 ]'

RC=0
"${COLLECT}" 0 a b 2> "${TMP_DIR}/usage.txt" || RC=$?
check "wrong argument count: exit 2 with a usage message" '[ "${RC}" -eq 2 ] && grep -q "^usage: semibin_collect_bins.sh" "${TMP_DIR}/usage.txt"'

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
