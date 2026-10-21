#!/bin/bash
#
# clean_read_count.sh <log> <samtools|bowtie2-paired|bowtie2-unpaired>
#
# Prints the number of clean reads a decontamination step wrote, read from
# the aligner's own stderr summary instead of a second pass over the gzipped
# output (performance review E8):
#
#   samtools          samtools fastq: "processed N reads" counts every record
#                     that passed -f/-F, "discarded M singletons" the paired
#                     records whose mate did not (sent to -s or dropped), so
#                     N - M reads were written to -1/-2 (or -0 for unpaired
#                     input, where M is 0)
#   bowtie2-paired    bowtie2 -1/-2 --un-conc: the pairs that "aligned
#                     concordantly 0 times", times two (reads in R1 + R2)
#   bowtie2-unpaired  bowtie2 -U --un: the reads that "aligned 0 times"
#
# Exits 1 when the expected line is missing, so a changed message format
# stops the task instead of reporting a wrong count.
#
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: clean_read_count.sh <log> <samtools|bowtie2-paired|bowtie2-unpaired>" >&2
    exit 2
fi

log="$1"
mode="$2"

first_match() {
    sed -nE "s/$1/\\1/p" "${log}" | head -n 1
}

case "${mode}" in
    samtools)
        processed=$(first_match '^\[M::bam2fq_mainloop\] processed ([0-9]+) reads$')
        singletons=$(first_match '^\[M::bam2fq_mainloop\] discarded ([0-9]+) singletons$')
        if [ -z "${processed}" ] || [ -z "${singletons}" ]; then
            echo "ERROR: no samtools fastq 'processed N reads' / 'discarded N singletons' lines in ${log}" >&2
            exit 1
        fi
        echo $(( processed - singletons ))
        ;;
    bowtie2-paired)
        pairs=$(first_match '^ *([0-9]+) \([0-9.]+%\) aligned concordantly 0 times$')
        if [ -z "${pairs}" ]; then
            echo "ERROR: no bowtie2 'aligned concordantly 0 times' line in ${log}" >&2
            exit 1
        fi
        echo $(( pairs * 2 ))
        ;;
    bowtie2-unpaired)
        reads=$(first_match '^ *([0-9]+) \([0-9.]+%\) aligned 0 times$')
        if [ -z "${reads}" ]; then
            echo "ERROR: no bowtie2 'aligned 0 times' line in ${log}" >&2
            exit 1
        fi
        echo "${reads}"
        ;;
    *)
        echo "ERROR: unknown mode '${mode}'" >&2
        exit 2
        ;;
esac
