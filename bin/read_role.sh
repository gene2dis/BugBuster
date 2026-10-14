#!/usr/bin/env bash
#
# read_role.sh <fastq>
#
# Prints the read role of a FASTQ file from its name: R1, R2 or S (singleton /
# unpaired), or nothing when the name carries none. Co-assembly pools every
# sample's reads into one list, so MEGAHIT and BOWTIE2_SAMTOOLS sort them by
# name (BWAMEM3_SAMTOOLS applies the same rule in Groovy).
#
# The extension is stripped and the name split on '_' and '.'; the LAST token
# that names a role wins (R1 or 1, R2 or 2, singleton or unpaired, any case).
# A sample id can therefore hold role-like parts: S_1.2_R2_clean -> R2 (a
# substring match on "_1." put it in R1), unpaired_s2_R1_clean -> R1,
# S_L001_R1_001 -> R1 ("001" is not a role).
set -euo pipefail

name=$(basename "$1")
for ext in .gz .bz2 .fastq .fq; do
    name=${name%"$ext"}
done

IFS='_.' read -r -a tokens <<< "$name"
for (( i = ${#tokens[@]} - 1; i >= 0; i-- )); do
    case "${tokens[i],,}" in
        r1|1) echo R1; exit 0 ;;
        r2|2) echo R2; exit 0 ;;
        singleton|unpaired) echo S; exit 0 ;;
    esac
done
