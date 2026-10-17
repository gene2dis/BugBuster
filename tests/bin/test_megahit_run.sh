#!/bin/bash
#
# Regression tests for bin/megahit_run.sh (MEGAHIT). Without -m, MEGAHIT
# sized its graph build to 90% of the host's RAM (1.9 TB on a 2 TB server)
# while docker capped the task at 72 GB; the kernel killed megahit_core, the
# wrapper exited 247, which is not in the retry list, and the co-assembly
# failed without a retry. Checks, with a stand-in megahit on PATH:
#   - -m is 90% of the given task memory and comes before the MEGAHIT
#     arguments (a -m in ext.args wins); '-' passes no -m;
#   - exit 247, and exit 1 with megahit_core's "No enough memory" in
#     <out dir>/log, become 137; other failures keep their status; success
#     stays 0;
#   - the out dir is found from -o and --out-dir;
#   - bad arguments are refused;
# and, with docker, that the real MEGAHIT (pinned image) runs through it and
# logs the memory budget.
#
# Requirements: bash; docker for the last check (skipped without it).
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
RUN="${REPO_DIR}/bin/megahit_run.sh"
MEGAHIT_IMG="quay.io/biocontainers/megahit:1.2.9--h43eeafb_4"

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

echo "=== megahit_run.sh test suite ==="

# Stand-in megahit: records its arguments, writes FAKE_LOG to <out dir>/log
# and exits FAKE_RC
mkdir -p "${TMP_DIR}/bin"
cat > "${TMP_DIR}/bin/megahit" <<'SH'
#!/bin/bash
printf '%s\n' "$@" > "${FAKE_ARGS}"
out=""
while [ $# -gt 0 ]; do
    case "$1" in -o|--out-dir) out="$2"; shift ;; --out-dir=*) out="${1#--out-dir=}" ;; esac
    shift
done
if [ -n "${out}" ] && [ -n "${FAKE_LOG:-}" ]; then
    mkdir -p "${out}" && echo "${FAKE_LOG}" > "${out}/log"
fi
exit "${FAKE_RC:-0}"
SH
chmod +x "${TMP_DIR}/bin/megahit"
export FAKE_ARGS="${TMP_DIR}/args.txt"

# fake <rc> <log text> <megahit_run.sh arguments...>: prints the exit status
fake() {
    local rc="$1" log="$2"
    shift 2
    rm -rf "${TMP_DIR}/w" && mkdir "${TMP_DIR}/w"
    ( cd "${TMP_DIR}/w" && PATH="${TMP_DIR}/bin:${PATH}" FAKE_RC="${rc}" FAKE_LOG="${log}" \
        "${RUN}" "$@" 2>/dev/null ) && echo 0 || echo $?
}

args_are() { [ "$(tr '\n' ' ' < "${FAKE_ARGS}")" = "$1" ]; }

# Memory budget
status=$(fake 0 "" 1000 -1 a_R1.fq -2 a_R2.fq -t 4 -o out)
check "success exits 0" [ "${status}" = 0 ]
check "-m is 90% of the task memory, before the MEGAHIT arguments" \
    args_are "-m 900 -1 a_R1.fq -2 a_R2.fq -t 4 -o out "
fake 0 "" 77309411328 -o out > /dev/null
check "-m for 72 GB is 69578470195 bytes" args_are "-m 69578470195 -o out "
fake 0 "" - -1 a_R1.fq -o out > /dev/null
check "'-' passes no -m" args_are "-1 a_R1.fq -o out "
fake 0 "" 1000 -m 0.5 -o out > /dev/null
check "a -m in the MEGAHIT arguments comes after the wrapper's (it wins)" args_are "-m 900 -m 0.5 -o out "

# Out of memory -> 137
check "exit 247 (megahit_core killed) becomes 137" [ "$(fake 247 "" 1000 -o out)" = 137 ]
check "exit 1 with 'No enough memory' in <out>/log becomes 137" \
    [ "$(fake 1 "FATAL sorting/base_engine.cpp: No enough memory" 1000 -o out)" = 137 ]
check "the log is found with --out-dir" \
    [ "$(fake 1 "No enough memory" 1000 --out-dir out)" = 137 ]
check "the log is found with --out-dir=" \
    [ "$(fake 1 "No enough memory" 1000 --out-dir=out)" = 137 ]
check "other failures keep their status (1)" [ "$(fake 1 "some other error" 1000 -o out)" = 1 ]
check "other failures keep their status (2)" [ "$(fake 2 "" 1000 -o out)" = 2 ]

# Bad arguments
check "a memory that is not a number of bytes is refused" [ "$(fake 0 "" 72GB -o out)" = 2 ]
check "too few arguments are refused" [ "$(fake 0 "" 1000)" = 2 ]

# Real MEGAHIT in the pinned image
real_megahit() {
    mkdir -p "${TMP_DIR}/real"
    cp "${REPO_DIR}/tests/data/reads/test_sample1_R1.fastq.gz" "${REPO_DIR}/tests/data/reads/test_sample1_R2.fastq.gz" "${TMP_DIR}/real/"
    docker run --rm -u "$(id -u):$(id -g)" \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" -v "${TMP_DIR}/real:${TMP_DIR}/real" -w "${TMP_DIR}/real" \
        "${MEGAHIT_IMG}" /pipeline_bin/megahit_run.sh 2000000000 \
        -1 test_sample1_R1.fastq.gz -2 test_sample1_R2.fastq.gz -t 2 -o real_assembly > /dev/null 2>&1 &&
    [ -f "${TMP_DIR}/real/real_assembly/final.contigs.fa" ] &&
    grep -q "Memory used: 1800000000" "${TMP_DIR}/real/real_assembly/log"
}
if command -v docker > /dev/null 2>&1; then
    check "real MEGAHIT runs through it with a 90% budget (Memory used: 1800000000)" real_megahit
else
    echo "- docker not found: real MEGAHIT check skipped"
fi

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
