#!/bin/bash
#
# Download and unpack a Kraken2 database tarball.
# Usage: kraken_db_reformat.sh <database-url>
#
set -euo pipefail

wget --tries=5 --continue --waitretry=10 -O k2_db.tar.gz \
    "${1:?usage: kraken_db_reformat.sh <database-url>}"

mkdir k2_ref_db
tar -xzf k2_db.tar.gz -C k2_ref_db

# A usable Kraken2 database always contains these three files
for f in hash.k2d opts.k2d taxo.k2d; do
    if [ ! -f "k2_ref_db/${f}" ]; then
        echo "ERROR: k2_ref_db/${f} missing after extraction" >&2
        exit 1
    fi
done
rm -f k2_db.tar.gz
