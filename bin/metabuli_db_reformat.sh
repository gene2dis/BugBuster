#!/bin/bash
#
# Download and unpack a Metabuli database (refseq_standard: an 80 GB download
# that unpacks to 110 GB; the tarball is deleted after extraction).
# Usage: metabuli_db_reformat.sh <database-tarball-url>
#
set -euo pipefail

wget --tries=5 --continue --waitretry=10 -O metabuli_db.tar.gz \
    "${1:?usage: metabuli_db_reformat.sh <database-tarball-url>}"

mkdir metabuli_extract
tar -xzf metabuli_db.tar.gz -C metabuli_extract
rm -f metabuli_db.tar.gz

# The tarball holds one database directory (or the files at its top level):
# the directory with db.parameters is the database
params_file=$(find metabuli_extract -maxdepth 3 -type f -name db.parameters | head -n 1)
if [ -z "${params_file}" ]; then
    echo "ERROR: no db.parameters found after extracting the Metabuli database" >&2
    exit 1
fi
mv "$(dirname "${params_file}")" metabuli_db
rm -rf metabuli_extract

for f in diffIdx info taxonomyDB; do
    if [ ! -s "metabuli_db/${f}" ]; then
        echo "ERROR: Metabuli database file ${f} missing or empty" >&2
        exit 1
    fi
done
# BLOBTOOLS reads the taxdump the database was built with
for f in nodes.dmp names.dmp; do
    if [ -z "$(find metabuli_db -maxdepth 2 -type f -name "${f}" | head -n 1)" ]; then
        echo "ERROR: the Metabuli database has no ${f} (its taxdump is needed by BlobTools)" >&2
        exit 1
    fi
done
