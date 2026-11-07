#!/bin/bash
#
# Download, verify and index a pinned CheckV database release (v1.5: a 1.68 GB
# download, 6.4 GB unpacked with its DIAMOND index). CheckV's own
# download_database takes whatever CURRENT_RELEASE.txt names, so the registry
# pins the file instead. CheckV publishes no checksum: the download's byte
# count, the release named in its README.txt and the files CheckV checks
# before a run are verified instead. The DIAMOND index is built as
# download_database builds it, with the image's DIAMOND.
# Usage: checkv_db_reformat.sh <database-tarball-url> <bytes> <version>
#
set -euo pipefail

usage="usage: checkv_db_reformat.sh <database-tarball-url> <bytes> <version>"
db_url="${1:?${usage}}"
db_bytes="${2:?${usage}}"
db_version="${3:?${usage}}"

wget --tries=5 --continue --waitretry=10 --no-verbose -O checkv_db.tar.gz "${db_url}"
found_bytes=$(stat -c %s checkv_db.tar.gz)
if [ "${found_bytes}" != "${db_bytes}" ]; then
    echo "ERROR: CheckV database download is ${found_bytes} bytes, expected ${db_bytes}" >&2
    exit 1
fi

mkdir checkv_extract
tar -xzf checkv_db.tar.gz -C checkv_extract
rm -f checkv_db.tar.gz
# The release tarball holds one checkv-db-v<version>/ directory
readme=$(find checkv_extract -maxdepth 2 -type f -name README.txt | head -n 1)
if [ -z "${readme}" ]; then
    echo "ERROR: no README.txt found after extracting the CheckV database" >&2
    exit 1
fi
mv "$(dirname "${readme}")" checkv_db
rm -rf checkv_extract

# README.txt opens with the newest release: '## Verion 1.5 (Jan 10, 2023)' (sic)
found_version=$(sed -n 's/^## Ver[a-z]*on \([0-9.]*\).*/\1/p' checkv_db/README.txt | head -n 1)
if [ "${found_version}" != "${db_version}" ]; then
    echo "ERROR: CheckV database version ${found_version}, expected ${db_version}" >&2
    exit 1
fi

diamond makedb \
    --in checkv_db/genome_db/checkv_reps.faa \
    --db checkv_db/genome_db/checkv_reps \
    --threads 1 \
    > /dev/null 2> checkv_db/genome_db/checkv_reps.log

for f in genome_db/checkv_reps.dmnd genome_db/checkv_reps.tsv genome_db/checkv_error.tsv hmm_db/checkv_hmms.tsv; do
    if [ ! -s "checkv_db/${f}" ]; then
        echo "ERROR: CheckV database file ${f} missing or empty" >&2
        exit 1
    fi
done
if [ ! -d checkv_db/hmm_db/checkv_hmms ]; then
    echo "ERROR: CheckV database directory hmm_db/checkv_hmms missing" >&2
    exit 1
fi
