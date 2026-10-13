#!/bin/bash
#
# Filter bwa-mem3 read -> contig alignments by identity and aligned length.
# Usage: bwa-mem3 mem ... | filter_alignments.sh <min_identity> <min_aligned_length> | samtools fixmate ...
#
# bwa mem aligns locally, so reads from genomes that are absent from (or only
# related to) the reference still place short or divergent clipped hits. In
# the bowtie2 vs bwa-mem3 benchmark (docs/benchmarks/bwamem3_vs_bowtie2.md)
# these inflated gene counts and bin depth in samples the assembly covers
# poorly; a MAPQ cut does not remove them, because such a hit is often the
# read's only placement. bowtie2 end-to-end rejected most of them.
#
# Reads SAM on stdin, writes SAM on stdout:
#   - header lines pass through;
#   - secondary and supplementary records are dropped (one record per read);
#   - unmapped records pass through;
#   - a mapped record is kept when it has >= <min_aligned_length> aligned
#     bases (CIGAR M/=/X) and an identity >= <min_identity>, with
#     identity = 1 - NM / (M + I + D);
#   - a failing record is turned into an unmapped one (flag 0x4 set, proper
#     pair 0x2 cleared, MAPQ 0, CIGAR *, TLEN 0, optional fields dropped), so
#     read totals and mapping rates still add up. RNAME/POS are kept; pipe the
#     output through `samtools fixmate` to update the mate fields.
# Both thresholds 0 keep every primary alignment.
#
# The input must keep a read's records together (bwa mem output does), which
# `samtools fixmate` also needs. Runs with POSIX awk (mawk in the container).
#
set -euo pipefail

min_identity="${1:?usage: filter_alignments.sh <min_identity> <min_aligned_length>}"
min_aligned_length="${2:?usage: filter_alignments.sh <min_identity> <min_aligned_length>}"

awk -F'\t' -v min_id="${min_identity}" -v min_len="${min_aligned_length}" '
    BEGIN {
        OFS = "\t"
        min_id += 0
        min_len += 0
        if (min_id < 0 || min_id > 1 || min_len < 0) {
            print "filter_alignments.sh: min_identity must be in [0, 1] and min_aligned_length >= 0" > "/dev/stderr"
            bad_args = 1
            exit 1
        }
    }
    /^@/ { print; next }
    # Bit tests by arithmetic: mawk has no and()/or()
    function has(flag, bit) { return int(flag / bit) % 2 }
    {
        flag = $2 + 0
        if (has(flag, 256) || has(flag, 2048)) next   # secondary, supplementary
        if (has(flag, 4)) { print; next }             # unmapped

        m = 0; ins = 0; del = 0
        cigar = $6
        while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
            n = substr(cigar, 1, RLENGTH - 1) + 0
            op = substr(cigar, RLENGTH, 1)
            if (op == "M" || op == "=" || op == "X") m += n
            else if (op == "I") ins += n
            else if (op == "D") del += n
            cigar = substr(cigar, RLENGTH + 1)
        }
        nm = 0
        for (i = 12; i <= NF; i++) {
            if (substr($i, 1, 5) == "NM:i:") { nm = substr($i, 6) + 0; break }
        }
        span = m + ins + del
        identity = span > 0 ? 1 - nm / span : 0
        if (m >= min_len && identity >= min_id) { print; next }

        flag += 4
        if (has(flag, 2)) flag -= 2
        print $1, flag, $3, $4, 0, "*", $7, $8, 0, $10, $11
    }
    END { if (bad_args) exit 1 }
'
