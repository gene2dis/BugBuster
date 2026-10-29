#!/bin/bash
#
# Download, verify and unpack a pinned geNomad database release (v1.9: a
# 0.84 GB download, 2.2 GB unpacked). geNomad's own download-database always
# fetches the latest release, so the registry pins the Zenodo file instead.
# Usage: genomad_db_reformat.sh <database-tarball-url> <md5> <version>
#
set -euo pipefail

db_url="${1:?usage: genomad_db_reformat.sh <database-tarball-url> <md5> <version>}"
db_md5="${2:?usage: genomad_db_reformat.sh <database-tarball-url> <md5> <version>}"
db_version="${3:?usage: genomad_db_reformat.sh <database-tarball-url> <md5> <version>}"

wget --tries=5 --continue --waitretry=10 --no-verbose -O genomad_db.tar.gz "${db_url}"
echo "${db_md5}  genomad_db.tar.gz" | md5sum -c -

mkdir genomad_extract
tar -xzf genomad_db.tar.gz -C genomad_extract
rm -f genomad_db.tar.gz

# The release tarball holds one genomad_db/ directory
version_file=$(find genomad_extract -maxdepth 2 -type f -name version.txt | head -n 1)
if [ -z "${version_file}" ]; then
    echo "ERROR: no version.txt found after extracting the geNomad database" >&2
    exit 1
fi
mv "$(dirname "${version_file}")" genomad_db
rm -rf genomad_extract

found_version=$(head -n 1 genomad_db/version.txt)
if [ "${found_version}" != "${db_version}" ]; then
    echo "ERROR: geNomad database version ${found_version}, expected ${db_version}" >&2
    exit 1
fi
for f in genomad_db genomad_db.dbtype genomad_marker_metadata.tsv; do
    if [ ! -s "genomad_db/${f}" ]; then
        echo "ERROR: geNomad database file ${f} missing or empty" >&2
        exit 1
    fi
done
