#!/bin/bash
#
# Compress the sourmash lineage CSVs staged into the working directory.
#
set -euo pipefail

# No *.csv staged makes gzip fail on the literal glob — that is a real input error
for f in *.csv; do
    gzip "${f}"
done
