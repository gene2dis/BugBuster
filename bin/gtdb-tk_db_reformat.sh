#!/bin/bash
#
# Download and unpack the GTDB-Tk reference data package into gtdbtk_db/
# (the package's release directory, e.g. release232/, renamed so the
# process output and its storeDir check have a fixed name).
# Usage: gtdb-tk_db_reformat.sh <package-url>
#
set -euo pipefail

wget --tries=5 --continue --waitretry=10 -O gtdbtk_data.tar.gz \
    "${1:?usage: gtdb-tk_db_reformat.sh <package-url>}"

mkdir gtdbtk_db
tar -xzf gtdbtk_data.tar.gz -C gtdbtk_db --strip-components=1

# The package extracts into one release directory; verify it is non-empty
n_files=$(find gtdbtk_db -type f | wc -l)
if [ "${n_files}" -eq 0 ]; then
    echo "ERROR: no files extracted from gtdbtk_data.tar.gz" >&2
    exit 1
fi
rm -f gtdbtk_data.tar.gz
