#!/bin/bash
#
# Regenerates the Nonpareil tables of tests/data/nonpareil/ (tests/bin/test_results_summary.sh) from the .npo
# curves committed there, with bin/nonpareil_fit.R in the pinned nonpareil image (the fit is deterministic):
#   sampleA.npo  Nonpareil 3.5.5 k-mer curve of R1 of a real 2 M-pair sample's clean reads (a model)
#   sampleB.npo  the same of tests/data/reads/test_sample1_R1.fastq.gz (500 reads: no model, Nd and LR* blank)
# Each gives <sample>_nonpareil.tsv and <sample>_nonpareil_curve.tsv next to it.
#
set -euo pipefail
REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
IMG=$(grep -o "'quay.io/biocontainers/nonpareil:[^']*'" "${REPO_DIR}/modules/local/nonpareil/main.nf" | tr -d "'")
docker run --rm -u "$(id -u):$(id -g)" -v "${REPO_DIR}:${REPO_DIR}" -w "${REPO_DIR}/tests/data/nonpareil" "${IMG}" sh -c "
    for npo in *.npo; do
        Rscript ${REPO_DIR}/bin/nonpareil_fit.R \"\${npo}\" \"\${npo%.npo}\"
    done"
