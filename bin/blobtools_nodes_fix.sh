#!/bin/sh
#
# nodes.dmp for BlobTools 1.1.1, which looks the top rank up as
# 'superkingdom'. NCBI renamed that rank 'domain' in 2025 and made Viruses
# (10239) an 'acellular root' with realms below it, so with a current taxdump
# BlobTools wrote 'undef' as every contig's superkingdom. Both get
# 'superkingdom' back (not realm: that would put realm names such as
# Riboviria in the column); a taxdump that already uses 'superkingdom' is
# unchanged. Writes the fixed nodes.dmp to stdout.
# Usage: blobtools_nodes_fix.sh <nodes.dmp>
#
set -eu

# The rank is the third field; in a three-column file it still carries the
# line's closing "\t|"
awk 'BEGIN { FS = "\t[|]\t"; OFS = "\t|\t" }
     ($3 ~ /^domain(\t[|])?$/ || $1 == "10239") { sub(/^[^\t]*/, "superkingdom", $3) }
     { print }' "${1:?usage: blobtools_nodes_fix.sh <nodes.dmp>}"
