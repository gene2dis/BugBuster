#!/bin/bash
#
# Regression tests for bin/read_role.sh, the read-role rule MEGAHIT and
# BOWTIE2_SAMTOOLS use to split pooled co-assembly reads into R1 / R2 /
# singleton lists (BWAMEM3_SAMTOOLS applies the same rule in Groovy, covered
# by tests/modules/bwamem3_samtools.nf.test). A substring match on "_1."
# put the R2 file of a sample id such as S_1.2 into R1.
#
# Requirements: bash only.
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
ROLE="${REPO_DIR}/bin/read_role.sh"
PASS=0
FAIL=0

expect() {
    local name="$1" want="$2" got
    got=$("${ROLE}" "${name}")
    if [ "${got}" = "${want}" ]; then
        echo "✓ ${name} -> ${want:-<none>}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${name}: expected '${want:-<none>}', got '${got:-<none>}'"
        FAIL=$((FAIL + 1))
    fi
}

# Post-QC names (<sample>_{R1,R2,Singleton}_clean.fastq.gz)
expect "s1_R1_clean.fastq.gz"                R1
expect "s1_R2_clean.fastq.gz"                R2
expect "s1_Singleton_clean.fastq.gz"         S
# Sample ids holding role-like parts: the last role token wins
expect "S_1.2_R1_clean.fastq.gz"             R1
expect "S_1.2_R2_clean.fastq.gz"             R2
expect "S_2.1_R1_clean.fastq.gz"             R1
expect "unpaired_s2_R1_clean.fastq.gz"       R1
expect "S_1.2_Singleton_clean.fastq.gz"      S
expect "/work/ab/cd/sim.ARG_R2_clean.fastq.gz" R2
# Raw input names (--quality_control false)
expect "sample_1.fastq.gz"                   R1
expect "sample_2.fq.gz"                      R2
expect "S_1.2_1.fastq.gz"                    R1
expect "S_1.2_2.fastq"                       R2
expect "S_L001_R1_001.fastq.gz"              R1
expect "S_L001_R2_001.fastq.gz"              R2
expect "sample.R2.fastq.gz"                  R2
expect "sample_unpaired.fastq.gz"            S
# No role in the name
expect "sample.fastq.gz"                     ""

echo
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
