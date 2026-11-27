#!/usr/bin/env bash
#
# semibin_collect_bins.sh <status> <bins_dir> <out_dir> <log> <sample>
#
# After a SemiBin2 binning command (SEMIBIN, both script branches): move its
# bin directory to the module's output, or write an empty output with a note
# when SemiBin2 found no bins, so Binette and the rest of the sample go on.
#   status    SemiBin2's exit status
#   bins_dir  the bin directory SemiBin2 writes (output_bins under its -o)
#   out_dir   the module's output directory (<sample>_semibin_output_bins)
#   log       SemiBin2's captured stdout + stderr
#   sample    the sample (or co-assembly) id, for the warning
#
# SemiBin2 fails when it finds no bins: the marker-gene search of its
# reclustering step gets an empty ORF file ("Number of bins prior to
# reclustering: 0"). SemiBin2 2.6.0 also stops, at its first step, when no
# contig reaches 2,500 bp ("but all are shorter than 2500 basepairs") or
# 4,000 bp, because it cannot build must-link pairs ("no must-link pairs can
# be generated"). Any other failure exits with SemiBin2's status. The three
# messages are SemiBin2 2.6.0's: check them when the image changes
# (tests/bin/test_semibin_collect_bins.sh). One copy for both modes (issue #30).
#
set -euo pipefail

if [ "$#" -ne 5 ]; then
    echo "usage: semibin_collect_bins.sh <status> <bins_dir> <out_dir> <log> <sample>" >&2
    exit 2
fi
status=$1
bins_dir=$2
out_dir=$3
log=$4
sample=$5

if [ "${status}" -eq 0 ]; then
    mv "${bins_dir}" "${out_dir}"
elif grep -q -e 'Number of bins prior to reclustering: 0' -e 'no must-link pairs can be generated' \
        -e 'but all are shorter than 2500 basepairs' "${log}"; then
    echo "WARNING: SemiBin found no bins for ${sample}." >&2
    mkdir -p "${out_dir}"
    echo "SemiBin found no bins" > "${out_dir}/SKIPPED.txt"
else
    exit "${status}"
fi
