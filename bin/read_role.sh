#!/usr/bin/env bash
#
# read_role.sh <fastq>
#
# Prints the read role of a FASTQ file from its name: R1, R2 or S (singleton /
# unpaired), or nothing when the name carries none. Co-assembly pools every
# sample's reads into one list, and MEGAHIT, BOWTIE2_SAMTOOLS and
# BWAMEM3_SAMTOOLS all sort them by name with this script.
#
# The extension is stripped and the name split on '_' and '.'. The LAST
# explicit role token wins (R1, R2, singleton, unpaired; any case); only a
# name without one falls back to its last bare 1 or 2. So a sample id can
# hold role-like parts, and chunk numbers do not override the role:
#   S_1.2_R2_clean -> R2   (a substring match on "_1." put it in R1)
#   unpaired_s2_R1_clean -> R1
#   S_L001_R1_001 -> R1    ("001" is not a role)
#   S1_R1_2 -> R1          (the trailing 2 is a chunk number)
#   S_1.2_1 -> R1          (no explicit token: last bare 1/2)
set -euo pipefail

name=$(basename "$1")
for ext in .gz .bz2 .fastq .fq; do
    name=${name%"$ext"}
done

IFS='_.' read -r -a tokens <<< "$name"
for (( i = ${#tokens[@]} - 1; i >= 0; i-- )); do
    case "${tokens[i],,}" in
        r1) echo R1; exit 0 ;;
        r2) echo R2; exit 0 ;;
        singleton|unpaired) echo S; exit 0 ;;
    esac
done
for (( i = ${#tokens[@]} - 1; i >= 0; i-- )); do
    case "${tokens[i]}" in
        1) echo R1; exit 0 ;;
        2) echo R2; exit 0 ;;
    esac
done
