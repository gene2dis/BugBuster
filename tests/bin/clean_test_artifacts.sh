#!/bin/bash
#
# Remove the local artifacts that nextflow / nf-test runs leave in the repo
# root. All of these are gitignored, so they never enter the repository —
# this script just reclaims the disk and keeps the working copy tidy.
#
# Deliberately NOT removed (real run outputs, also gitignored):
#   results/  results_stub/  databases/
#
# Usage: tests/bin/clean_test_artifacts.sh
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"

cd "${REPO_DIR}"

removed=0
for target in work .nextflow .nf-test test_results null; do
    if [ -e "${target}" ]; then
        rm -rf "${target}"
        echo "removed ${target}/"
        removed=$((removed + 1))
    fi
done

for pattern in '.nextflow.log' '.nextflow.log.*' '.nf-test.log' '.nf-test.log.*'; do
    for f in ${pattern}; do
        if [ -f "${f}" ]; then
            rm -f "${f}"
            echo "removed ${f}"
            removed=$((removed + 1))
        fi
    done
done

if [ "${removed}" -eq 0 ]; then
    echo "nothing to clean"
else
    echo "cleaned ${removed} test artifact(s) from ${REPO_DIR}"
fi
