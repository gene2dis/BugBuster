#!/bin/bash
#
# Download and unpack the CheckM2 DIAMOND database.
# Usage: checkm2_db_reformat.sh <database-url>
#
set -euo pipefail

wget --tries=5 --continue --waitretry=10 -O checkm2_database.tar.gz \
    "${1:?usage: checkm2_db_reformat.sh <database-url>}"

tar -xvzf checkm2_database.tar.gz
find . -mindepth 2 -type f -name '*.dmnd' | xargs -I {} mv {} .

# Verify the database landed before cleaning up, so exit status reflects the DB
if ! ls ./*.dmnd > /dev/null 2>&1; then
    echo "ERROR: no .dmnd file found after extracting checkm2_database.tar.gz" >&2
    exit 1
fi

rm -f checkm2_database.tar.gz CONTENTS.json
rm -rf CheckM2_database
