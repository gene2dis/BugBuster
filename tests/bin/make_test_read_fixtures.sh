#!/bin/bash
#
# Regenerate tests/data/reads/test_sample{1,2}_R{1,2}.fastq.gz: the first 500
# read pairs of the nf-core test-datasets Bacteroides fragilis reads (MIT
# licence) the test profile uses (assets/test_samplesheet.csv). The stub and
# validation nf-tests read these instead of downloading from GitHub, whose
# raw-file downloads failed intermittently on CI runners (2026-10-09). Read
# content does not matter to a stub run; a real run still uses the full files.
#
# Output is byte-stable (gzip -n). Requires curl and network access.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUT_DIR="$( cd "${SCRIPT_DIR}/../data/reads" && pwd )"
BASE="https://github.com/nf-core/test-datasets/raw/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq"
N_PAIRS=500

for s in 1 2; do
    for r in 1 2; do
        out="${OUT_DIR}/test_sample${s}_R${r}.fastq.gz"
        # head closes the pipe early, so the upstream curl/gzip exit non-zero by
        # design; the line count below is the check
        set +o pipefail
        curl -sSfL "${BASE}/test${s}_${r}.fastq.gz" 2>/dev/null | gzip -cd 2>/dev/null \
            | head -n $((N_PAIRS * 4)) | gzip -n > "${out}"
        set -o pipefail
        lines=$(gzip -cd "${out}" | wc -l)
        if [ "${lines}" -ne $((N_PAIRS * 4)) ]; then
            echo "ERROR: ${out} has ${lines} lines, expected $((N_PAIRS * 4))" >&2
            exit 1
        fi
        echo "wrote ${out} (${N_PAIRS} reads)"
    done
done
