#!/bin/bash
#
# Every shell test suite runs in CI (roadmap Phase 25, RV27; review finding
# R8-15): ci.yml has one step per tests/bin/test_*.sh, written by hand, so a
# new suite used to sit unused until someone added its step. Fails when
#   - a tests/bin/test_*.sh is not called in .github/workflows/ci.yml,
#   - it is called where its failure cannot fail CI ("|| true" after it,
#     continue-on-error: true or if: false on its step or job), or
#   - ci.yml calls a tests/bin/test_*.sh that does not exist.
# A suite that runs elsewhere on purpose would be listed in NOT_IN_CI with
# the reason (none today).
#
# Requires python3 with PyYAML (CI: apt python3-yaml).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
CI="${CI_WORKFLOW:-${REPO_DIR}/.github/workflows/ci.yml}"
NOT_IN_CI=()

PASS=0
FAIL=0
ok() { echo "✓ $1"; PASS=$((PASS + 1)); }
ko() { echo "✗ $1"; FAIL=$((FAIL + 1)); }

echo "=== CI wiring of the tests/bin suites ==="

# The suites ci.yml runs, read from its steps (python3 with PyYAML, as on the
# CI runners): a call counts when a logical line of a step's run script starts
# with tests/bin/test_<name>.sh (after optional VAR=value assignments) and its
# failure fails the job: not followed by "||", and neither the step nor its
# job has continue-on-error: true or if: false
calls=$(python3 -I - "${CI}" <<'PY'
import re, sys
import yaml

wf = yaml.safe_load(open(sys.argv[1]))
off = lambda v: str(v).strip().lower() in ("false", "${{ false }}")
for job_id, job in (wf.get("jobs") or {}).items():
    job_quiet = job.get("continue-on-error") is True or off(job.get("if", ""))
    for step in job.get("steps") or []:
        quiet = job_quiet or step.get("continue-on-error") is True or off(step.get("if", ""))
        script = re.sub(r"\\\n", " ", step.get("run") or "")
        for line in script.splitlines():
            m = re.match(r"\s*(?:[A-Za-z_][A-Za-z0-9_]*=(?:\"[^\"]*\"|'[^']*'|\S*)\s+)*(tests/bin/test_[A-Za-z0-9_]+\.sh)\b(.*)", line)
            if m:
                state = "quiet" if quiet or "||" in m.group(2) else "ok"
                print(f"{m.group(1)}\t{state}\t{job_id}")
PY
) || { ko "could not read ${CI#"${REPO_DIR}/"} as YAML"; calls=""; }
called=$(printf '%s\n' "${calls}" | awk -F'\t' '$2 == "ok" {print $1}' | sort -u)
quieted=$(printf '%s\n' "${calls}" | awk -F'\t' '$2 == "quiet" {print $1}' | sort -u)
n_suites=0
for f in "${REPO_DIR}"/tests/bin/test_*.sh; do
    rel="tests/bin/$(basename "${f}")"
    n_suites=$((n_suites + 1))
    if printf '%s\n' "${NOT_IN_CI[@]+"${NOT_IN_CI[@]}"}" | grep -qx "${rel}"; then
        ok "${rel} (not in CI on purpose)"
    elif printf '%s\n' "${called}" | grep -qx "${rel}"; then
        ok "${rel} runs in ci.yml"
    elif printf '%s\n' "${quieted}" | grep -qx "${rel}"; then
        ko "${rel} runs in ${CI#"${REPO_DIR}/"} but cannot fail it (|| after the call, continue-on-error or if: false)"
    else
        ko "${rel} is not run by ${CI#"${REPO_DIR}/"}: add a step for it"
    fi
done
while read -r rel; do
    [ -z "${rel}" ] && continue
    if [ -f "${REPO_DIR}/${rel}" ]; then
        :
    else
        ko "ci.yml runs ${rel}, which does not exist"
    fi
done <<< "$(printf '%s\n%s\n' "${called}" "${quieted}")"
[ "${n_suites}" -gt 0 ] && ok "${n_suites} suites found" || ko "no tests/bin/test_*.sh found"

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
