#!/bin/bash
#
# Download and unpack the GTDB-Tk reference data package.
# Usage: gtdb-tk_db_reformat.sh <package-url>
#
set -euo pipefail

wget --tries=5 --continue --waitretry=10 -O gtdbtk_data.tar.gz \
    "${1:?usage: gtdb-tk_db_reformat.sh <package-url>}"

tar -xzf gtdbtk_data.tar.gz

# The package extracts into a release directory; verify it is non-empty
n_files=$(find . -mindepth 2 -type f | wc -l)
if [ "${n_files}" -eq 0 ]; then
    echo "ERROR: no files extracted from gtdbtk_data.tar.gz" >&2
    exit 1
fi
rm -f gtdbtk_data.tar.gz
