#!/bin/bash
#
# The settings each profile makes are in force in the resolved config (roadmap
# Phase 22, RV15; review R2-02, R2-03, R2-13), checked with `nextflow config
# -profile <p> -flat` on the Nextflow the job installs (CI runs it on the
# minimum version and the latest: the legacy and the strict config parser
# order blocks differently).
#
# Before Phase 22 the process {} rules and conf/base.config came after the
# profiles {} block: under the strict parser (Nextflow >= 26.04) they replaced
# slurm's executor and queue size (16 instead of 100) and the cloud profiles'
# retries. Placing them before the profiles in turn let the legacy parser's
# profile withLabel blocks replace the whole label block (cpus, memory, time),
# so the labels must merge with slurm's queue routing: checked below.
#
# Also checked, from the resolved config of every profile:
#   - every errorStrategy starts with the terminate-on-unschedulable prefix, and
#     one that retries listed exit codes also retries the cloud preemption codes
#     (a withLabel / withName rule replaces the default under any profile);
#   - every database process stored with storeDir has a publishDir of its own,
#     so conf/base.config's default rule never copies a database into the
#     results directory (RV14; review R2-01);
# and, from the files, every errorStrategy in conf/ (the institutional template
# is not a profile) carries the prefix.
#
# `nextflow config` takes no parameters, so --slurm_queue_size given on the
# command line is not checked here; it is read when the config is parsed, like
# --max_cpus, whose command-line path tests/bin/test_resource_limits.sh covers.
#
# Requirements: nextflow (plugins for the cloud profiles are downloaded on
# first use).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0
pass() { echo "✓ $1"; PASS=$((PASS + 1)); }
fail() { echo "✗ $1"; FAIL=$((FAIL + 1)); }

PREFIX="task.exitStatus == Integer.MAX_VALUE && task.executor == 'local' ? 'terminate' :"
PREEMPT="task.exitStatus == ['google-batch': 14, 'azurebatch': 50001][task.executor]"

resolve() {
    # resolve <profile|-> : the flat config, in ${TMP_DIR}/<name>.txt
    local name="${1//,/+}"
    local args=()
    [ "$1" != "-" ] && args=(-profile "$1")
    if ! (cd "${TMP_DIR}" && nextflow config "${args[@]}" -flat "${REPO_DIR}" > "${TMP_DIR}/${name}.txt" 2> "${TMP_DIR}/${name}.err"); then
        fail "nextflow config ${args[*]} failed: $(head -3 "${TMP_DIR}/${name}.err")"
        return 1
    fi
}

expect() {
    # expect <profile> <exact flat line>...
    local profile="$1"; shift
    local file="${TMP_DIR}/${profile//,/+}.txt" line ok=1
    for line in "$@"; do
        if ! grep -qxF -- "${line}" "${file}"; then
            fail "${profile}: missing '${line}'"
            ok=0
        fi
    done
    [ ${ok} -eq 1 ] && pass "${profile}: $# settings in force"
    return 0
}

ENGINES="conda docker singularity podman apptainer shifter charliecloud"
engine_on() {
    # engine_on <engine>: the lines of a profile that enables <engine> only
    local e lines=()
    for e in ${ENGINES}; do
        if [ "${e}" = "$1" ]; then lines+=("${e}.enabled = true"); else lines+=("${e}.enabled = false"); fi
    done
    printf '%s\n' "${lines[@]}"
}

PROFILES="- debug conda docker singularity podman apptainer local_docker slurm_singularity test aws gcp azure slurm low_disk slurm,singularity slurm,apptainer aws,docker gcp,docker azure,docker test,docker"

echo "=== Resolved profile config ($(nextflow -version 2>/dev/null | grep -o 'version [0-9.]*' | head -1)) ==="
for p in ${PROFILES}; do
    resolve "${p}" || true
done

# --- No profile: the local pool from the max_* defaults, no database storage
expect - "executor.name = 'local'" "executor.cpus = 16" "executor.memory = '128.GB'" \
    "executor.queueSize = 16" "params.databases_dir = null"

# --- Container engines (local_docker is the docker profile)
for p in docker local_docker; do
    mapfile -t lines < <(engine_on docker)
    expect "${p}" "${lines[@]}" "docker.runOptions = '-u \$(id -u):\$(id -g)'"
done
for e in singularity podman apptainer conda; do
    mapfile -t lines < <(engine_on "${e}")
    expect "${e}" "${lines[@]}"
done

# --- SLURM: executor and queue size of the profile, labels merged with its queues
SLURM_LINES=(
    "process.executor = 'slurm'"
    "executor.name = 'slurm'"
    "executor.queueSize = 100"
    "executor.submitRateLimit = '10 sec'"
    "process.maxRetries = 2"
    "process.queue = 'normal'"
    "process.'withLabel:process_single'.queue = 'normal'"
    "process.'withLabel:process_single'.time = { 4.h * task.attempt }"
    "process.'withLabel:process_high'.queue = 'normal'"
    "process.'withLabel:process_high'.cpus = 16"
    "process.'withLabel:process_high'.memory = { 72.GB * task.attempt }"
    "process.'withLabel:process_high'.time = { 72.h * task.attempt }"
)
for p in slurm slurm,singularity slurm,apptainer slurm_singularity; do
    expect "${p}" "${SLURM_LINES[@]}"
done
mapfile -t lines < <(engine_on singularity)
expect slurm_singularity "${lines[@]}" "params.config_profile_name = 'SLURM HPC'"

# --- Cloud executors: their executor and retries
for p in aws aws,docker; do
    expect "${p}" "process.executor = 'awsbatch'" "process.maxRetries = 3" "docker.enabled = true"
done
for p in gcp gcp,docker; do
    expect "${p}" "process.executor = 'google-batch'" "process.maxRetries = 3" "docker.enabled = true"
done
for p in azure azure,docker; do
    expect "${p}" "process.executor = 'azurebatch'" "process.maxRetries = 3" "docker.enabled = true"
done

# --- Test profile: its own pool, and databases_dir from its output
for p in test test,docker; do
    expect "${p}" "executor.cpus = 2" "executor.memory = '6.GB'" "executor.queueSize = 2" \
        "params.output = 'test_results'" "params.databases_dir = 'test_results/../databases'"
done

# --- Every profile: error strategies and database publishing
for p in ${PROFILES}; do
    file="${TMP_DIR}/${p//,/+}.txt"
    [ -s "${file}" ] || continue
    bad=$(grep -E "errorStrategy = " "${file}" | grep -vF -- "${PREFIX}" || true)
    [ -z "${bad}" ] || fail "${p}: errorStrategy without the terminate prefix: ${bad}"
    bad=$(grep -E "errorStrategy = .*130\.\.145" "${file}" | grep -vF -- "${PREEMPT}" || true)
    [ -z "${bad}" ] || fail "${p}: errorStrategy retrying listed codes without the preemption codes: ${bad}"
    stored=$(grep -oE "^process\.'withName:[^']+'\.storeDir" "${file}" | sed -E "s/\.storeDir$//" || true)
    missing=""
    for sel in ${stored}; do
        grep -qF "${sel}.publishDir" "${file}" || missing="${missing} ${sel}"
    done
    [ -z "${missing}" ] || fail "${p}: stored without a publishDir of its own (base.config's default would copy them):${missing}"
    [ -z "${bad}${missing}" ] && pass "${p}: error strategies and database publishing"
done

# --- conf/ files that are not profiles (templates included with -c)
for f in "${REPO_DIR}"/conf/*.config; do
    bad=$(grep -E "^\s*errorStrategy\s*=" "${f}" | grep -vF -- "${PREFIX}" || true)
    if [ -n "${bad}" ]; then fail "$(basename "${f}"): errorStrategy without the prefix: ${bad}"; fi
done
pass "conf/*.config: every errorStrategy carries the prefix"

echo ""
echo "Passed: ${PASS}, failed: ${FAIL}"
[ ${FAIL} -eq 0 ]
