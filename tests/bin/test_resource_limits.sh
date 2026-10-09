#!/bin/bash
#
# Regression test for design doc Q22 (T9 finding, 2026-10-09): max_* params set
# in a -c config file never reached process.resourceLimits, because
# nextflow.config copied the params into a plain map while it was parsed —
# after `-profile test` every task stayed clamped at 2 CPU / 6 GB however the
# -c file raised them. resourceLimits is now a closure read at task time; the
# local executor pool (conf/base.config) still copies max_* at parse time, so a
# -c file must set executor.cpus/memory too, and otherwise fails and exits: a
# local task that can never be scheduled ends the run ('terminate'; with the
# old 'finish' the run reported the error and then hung, Q22 reopened).
#
# Each case is a stub launch of the QC-only test profile under docker; the
# per-task limits are read from the `--cpu-shares` / `--memory` options
# Nextflow writes into each task's .command.run (1024 shares per CPU).
#   - -c (params + executor): tasks clamped at the raised 4 CPU / 12 GB
#   - -c (params only, raised above the pool): fails with "exceeds available" and
#     exits (a hang is caught by the 300 s timeout and fails the test)
#   - -params-file and --max_* on the command line: tasks at 4 CPU / 12 GB
#   - no override: tasks stay at the test profile's 2 CPU / 6 GB
#
# Requirements: nextflow and docker (as in the CI pipeline job).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"

TMP_DIR=$(mktemp -d)
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

echo "=== Resource limit override test suite ==="
echo ""

printf "params { max_cpus = 4; max_memory = '12.GB'; max_time = '2.h' }\nexecutor { cpus = 4; memory = '12.GB' }\n" \
    > "${TMP_DIR}/limits_with_executor.config"
printf "params { max_cpus = 4; max_memory = '12.GB'; max_time = '2.h' }\n" > "${TMP_DIR}/limits_params_only.config"
printf "max_cpus: 4\nmax_memory: '12.GB'\nmax_time: '2.h'\n" > "${TMP_DIR}/limits.yaml"
# Local reads (tests/bin/make_test_read_fixtures.sh): the stub launches need no download
READS="${REPO_DIR}/tests/data/reads"
printf "sample,r1,r2,s\ntest_sample1,%s,%s,\ntest_sample2,%s,%s,\n" \
    "${READS}/test_sample1_R1.fastq.gz" "${READS}/test_sample1_R2.fastq.gz" \
    "${READS}/test_sample2_R1.fastq.gz" "${READS}/test_sample2_R2.fastq.gz" > "${TMP_DIR}/samplesheet.csv"

PASS=0
FAIL=0
pass() { echo "✓ $1"; PASS=$((PASS + 1)); }
fail() { echo "✗ $1"; FAIL=$((FAIL + 1)); }

launch() {
    # launch <case-dir> <nextflow global opts> <run opts> — returns the run's exit status
    local dir="${TMP_DIR}/$1"
    mkdir -p "${dir}"
    # timeout: a run that reports an error but never exits must not pass (or
    # stall CI); exit status 124 = hung
    (cd "${dir}" && timeout 300 nextflow $2 run "${REPO_DIR}/main.nf" -profile test,docker -stub \
        --input "${TMP_DIR}/samplesheet.csv" --output results $3 > run.log 2>&1)
}

max_limits() {
    # max_limits <case-dir> — "<max cpu-shares> <max memory MB>" over every task
    local dir="${TMP_DIR}/$1"
    local shares mem
    shares=$(cat "${dir}"/work/*/*/.command.run | grep -o -- '--cpu-shares [0-9]*' | awk '{print $2}' | sort -n | tail -1)
    mem=$(cat "${dir}"/work/*/*/.command.run | grep -o -- '--memory [0-9]*m' | tr -dc '0-9\n' | sort -n | tail -1)
    echo "${shares:-none} ${mem:-none}"
}

expect_limits() {
    # expect_limits <desc> <case-dir> <cpu-shares> <memory MB>
    local got
    got=$(max_limits "$2")
    if [ "${got}" = "$3 $4" ]; then
        pass "$1 (max --cpu-shares $3, --memory $4m)"
    else
        fail "$1 — expected max '$3 $4', got '${got}'"
    fi
}

if launch c_executor "-c ${TMP_DIR}/limits_with_executor.config" ""; then
    pass "-c with params + executor: run succeeds"
else
    fail "-c with params + executor: run succeeds"; tail -5 "${TMP_DIR}/c_executor/run.log"
fi
expect_limits "-c with params + executor: tasks clamped at the raised limits" c_executor 4096 12288

rc=0
launch c_params_only "-c ${TMP_DIR}/limits_params_only.config" "" || rc=$?
if [ "${rc}" -eq 0 ]; then
    fail "-c with params only (above the 2-CPU test pool): expected a scheduling failure"
elif [ "${rc}" -eq 124 ]; then
    fail "-c with params only: reported the error but did not exit (hung, killed after 300 s)"
elif grep -q "exceeds available" "${TMP_DIR}/c_params_only/run.log"; then
    pass "-c with params only (above the 2-CPU test pool): fails loudly ('exceeds available')"
else
    fail "-c with params only: failed, but not with 'exceeds available'"; tail -5 "${TMP_DIR}/c_params_only/run.log"
fi

if launch params_file "" "-params-file ${TMP_DIR}/limits.yaml"; then
    pass "-params-file: run succeeds"
else
    fail "-params-file: run succeeds"; tail -5 "${TMP_DIR}/params_file/run.log"
fi
expect_limits "-params-file: tasks clamped at the raised limits" params_file 4096 12288

if launch cli "" "--max_cpus 4 --max_memory 12.GB --max_time 2.h"; then
    pass "--max_* on the command line: run succeeds"
else
    fail "--max_* on the command line: run succeeds"; tail -5 "${TMP_DIR}/cli/run.log"
fi
expect_limits "--max_* on the command line: tasks clamped at the raised limits" cli 4096 12288

if launch default "" ""; then
    pass "no override: run succeeds"
else
    fail "no override: run succeeds"; tail -5 "${TMP_DIR}/default/run.log"
fi
expect_limits "no override: tasks stay at the test profile's 2 CPU / 6 GB" default 2048 6144

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
