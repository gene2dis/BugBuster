#!/bin/bash
#
# Every shell test suite runs in CI (roadmap Phase 25, RV27; review finding
# R8-15): ci.yml has one step per tests/bin/test_*.sh, written by hand, so a
# new suite used to sit unused until someone added its step. Fails when
#   - a tests/bin/test_*.sh is not called in .github/workflows/ci.yml, or
#   - ci.yml calls a tests/bin/test_*.sh that does not exist.
# A suite that runs elsewhere on purpose would be listed in NOT_IN_CI with
# the reason (none today).
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

# The suites ci.yml runs: tests/bin/test_<name>.sh at the start of a command
called=$(grep -oE '^[[:space:]]*([A-Z_]+="[^"]*"[[:space:]]+)?tests/bin/test_[A-Za-z0-9_]+\.sh' "${CI}" \
    | grep -oE 'tests/bin/test_[A-Za-z0-9_]+\.sh' | sort -u)
n_suites=0
for f in "${REPO_DIR}"/tests/bin/test_*.sh; do
    rel="tests/bin/$(basename "${f}")"
    n_suites=$((n_suites + 1))
    if printf '%s\n' "${NOT_IN_CI[@]+"${NOT_IN_CI[@]}"}" | grep -qx "${rel}"; then
        ok "${rel} (not in CI on purpose)"
    elif printf '%s\n' "${called}" | grep -qx "${rel}"; then
        ok "${rel} runs in ci.yml"
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
done <<< "${called}"
[ "${n_suites}" -gt 0 ] && ok "${n_suites} suites found" || ko "no tests/bin/test_*.sh found"

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
