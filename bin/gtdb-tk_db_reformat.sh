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

mkdir gtdbtk_extract
tar -xzf gtdbtk_data.tar.gz -C gtdbtk_extract

# The data root is the directory holding markers/ (release232/, possibly
# behind a ./ prefix); GTDBTK_DATA_PATH points there. storeDir keeps a wrong
# layout for every later run, so check it here
markers=$(find gtdbtk_extract -maxdepth 3 -type d -name markers | head -n 1)
if [ -z "${markers}" ] || [ -z "$(find "${markers}" -type f | head -n 1)" ]; then
    echo "ERROR: no markers/ with files in gtdbtk_data.tar.gz (not a GTDB-Tk data package?)" >&2
    exit 1
fi
mv "$(dirname "${markers}")" gtdbtk_db
rm -rf gtdbtk_extract gtdbtk_data.tar.gz
