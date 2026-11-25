#!/bin/bash
#
# Convert a Prodigal/Pyrodigal GFF3 (optionally gzipped) to featureCounts SAF.
# Usage: gff2saf.sh <genes.gff[.gz]> <out.saf>
#
# Emits the SAF header (GeneID, Chr, Start, End, Strand) followed by one row
# per CDS feature. The GeneID is built as
# "<GFF column 1 seqid>_<genenum from the ID suffix>", the form of the FAA
# headers that eggNOG-mapper annotates, so featureCounts gene ids join directly
# against the annotation table (design doc Section 4.4). Pyrodigal (the
# pipeline's gene caller) already writes ID=<contig_name>_<genenum> in its GFF
# (e.g. ID=k141_0_5); Prodigal writes ID=<seqnum>_<genenum> (e.g. ID=1_5, the
# ordinal of the contig in the input FASTA), which is also the ID= inside
# Pyrodigal's FAA header comments. The rewrite gives the same GeneID for both
# (review R6-07).
#
# An empty gene set (e.g. an empty-contig sample) yields a header-only SAF.
# A CDS row without a parseable ID attribute is an error, not a silent skip.
#
set -euo pipefail

in_gff="${1:?usage: gff2saf.sh <genes.gff[.gz]> <out.saf>}"
out_saf="${2:?usage: gff2saf.sh <genes.gff[.gz]> <out.saf>}"

# gzip -cdf passes plain (non-gzipped) files through unchanged (repo idiom)
gzip -cdf "${in_gff}" | awk -F'\t' '
    BEGIN {
        OFS = "\t"
        print "GeneID", "Chr", "Start", "End", "Strand"
    }
    /^#/ { next }
    NF >= 9 && $3 == "CDS" {
        id = ""
        n = split($9, attrs, ";")
        for (i = 1; i <= n; i++) {
            if (attrs[i] ~ /^ID=/) {
                id = substr(attrs[i], 4)
            }
        }
        genenum = id
        sub(/^.*_/, "", genenum)
        if (id == "" || genenum !~ /^[0-9]+$/) {
            printf "ERROR: gff2saf.sh: missing or malformed ID attribute on GFF line %d: %s\n", NR, $9 > "/dev/stderr"
            exit 1
        }
        print $1 "_" genenum, $1, $4, $5, $7
    }
' > "${out_saf}"
