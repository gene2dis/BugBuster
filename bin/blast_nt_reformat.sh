#!/bin/bash
#
# Download and assemble the NCBI nt BLAST database.
# Usage: blast_nt_reformat.sh <ncbi-blast-db-base-url>   (e.g. https://ftp.ncbi.nlm.nih.gov/blast/db)
#
set -euo pipefail

base_url="${1:?usage: blast_nt_reformat.sh <ncbi-blast-db-base-url>}"

# NCBI publishes a JSON manifest of the volumes; use it instead of scraping the HTML listing
wget --tries=5 --waitretry=10 -qO nt-nucl-metadata.json "${base_url}/nt-nucl-metadata.json"
volumes=$(grep -o 'nt\.[0-9][0-9]*\.tar\.gz' nt-nucl-metadata.json | sort -u)
if [ -z "${volumes}" ]; then
    echo "ERROR: no nt volumes listed in ${base_url}/nt-nucl-metadata.json" >&2
    exit 1
fi
n_volumes=$(echo "${volumes}" | wc -l)
echo "Downloading ${n_volumes} nt volumes from ${base_url}"

# Download every volume plus its md5 file; a failure in any xargs job fails the script
echo "${volumes}" | xargs -n1 -P4 -I{} \
    wget --tries=5 --continue --waitretry=10 --no-verbose "${base_url}/{}" "${base_url}/{}.md5"

# Verify checksums before extracting anything
for vol in ${volumes}; do
    md5sum -c "${vol}.md5"
done

for vol in ${volumes}; do
    tar -xzf "${vol}"
done

# Verify the assembled database before moving it into place
for f in taxdb.btd taxdb.bti taxonomy4blast.sqlite3; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: expected file '${f}' missing after extraction" >&2
        exit 1
    fi
done
n_indices=$(find . -maxdepth 1 -name 'nt.*.nin' | wc -l)
if [ "${n_indices}" -ne "${n_volumes}" ]; then
    echo "ERROR: extracted ${n_indices} volume indices (nt.*.nin), expected ${n_volumes}" >&2
    exit 1
fi

rm -f nt.*.tar.gz nt.*.tar.gz.md5 nt-nucl-metadata.json
mkdir blast_nt_db
mv taxdb.btd taxdb.bti taxonomy4blast.sqlite3 blast_nt_db
mv nt.* blast_nt_db
