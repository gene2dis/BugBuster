#!/bin/bash
#
# Download the eggNOG 7 data files needed by eggNOG-mapper v3 (diamond mode).
# Usage: eggnog_db_reformat.sh <emapper-3.0-data-base-url>   (e.g. https://data.cgmlab.org/eggnog-mapper/emapper-3.0/data)
#
# Produces eggnog_db/ with the canonical emapper-3.0 data-dir layout:
#   eggnog.db, eggnog.db.fieldpresence.bin, eggnog.db.taxids.bin,
#   eggnog.taxa.db, eggnog.taxa.db.traverse.pkl, eggnog_proteins.dmnd,
#   go-basic.obo
# The files are served uncompressed (~44 GB total). The optional mmseqs/
# subdirectory is not fetched — the pipeline runs diamond mode only.
#
set -euo pipefail

base_url="${1:?usage: eggnog_db_reformat.sh <emapper-3.0-data-base-url>}"

files="eggnog.db eggnog.db.fieldpresence.bin eggnog.db.taxids.bin eggnog.taxa.db eggnog.taxa.db.traverse.pkl eggnog_proteins.dmnd go-basic.obo"

mkdir -p eggnog_db
for f in ${files}; do
    wget --tries=5 --continue --waitretry=10 --no-verbose -O "eggnog_db/${f}" "${base_url}/${f}"
done

# Verify the database landed before finishing, so exit status reflects the DB
for f in ${files}; do
    if [ ! -s "eggnog_db/${f}" ]; then
        echo "ERROR: expected file 'eggnog_db/${f}' missing or empty after download" >&2
        exit 1
    fi
done
