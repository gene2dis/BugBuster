#!/bin/bash
#
# Remove the local artifacts that nextflow / nf-test runs leave in the repo
# root. All of these are gitignored, so they never enter the repository —
# this script just reclaims the disk and keeps the working copy tidy.
#
# Removed: work/ .nextflow/ .nextflow.log* .nf-test/ .nf-test.log* test_results/,
# and databases/stub/: a -stub run of -profile test stores its placeholder
# decontamination index there (only -stub runs write under stub/).
#
# Deliberately NOT removed (real run outputs, also gitignored):
#   results/  results_stub/  the rest of databases/ (downloaded databases; the
#   full test run's databases/bwamem3_index/ shares its folder with a real run
#   that uses --custom_host_fasta, so it is left for you to delete)
#
# Usage: tests/bin/clean_test_artifacts.sh
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"

cd "${REPO_DIR}"

removed=0
for target in work .nextflow .nf-test test_results databases/stub; do
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
