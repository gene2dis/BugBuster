#!/bin/bash
#
# Download one half of a sylph reference: the genome database or the
# sylph-tax taxonomy metadata that matches it.
# Usage: sylph_db_reformat.sh db       <database-url>   (.syldb or .syl2db)
#        sylph_db_reformat.sh taxonomy <metadata-url>   (.tsv.gz)
#
# The file lands in sylph_db/ under its URL name, which carries the release
# (e.g. gtdb-r232-c200-dbv2.syl2db, gtdb_r232_metadata.tsv.gz) and is
# recorded as provenance by SYLPH_PROFILE / SYLPH_TAXPROF. The two halves
# run as separate tasks (FORMAT_SYLPH_DB / FORMAT_SYLPH_TAXONOMY) with their
# own checks and retries; --custom_sylph_db and --custom_sylph_taxonomy
# replace both together.
#
set -euo pipefail

usage="usage: sylph_db_reformat.sh <db|taxonomy> <url>"
kind="${1:?${usage}}"
url="${2:?${usage}}"
name="${url##*/}"

case "${kind}" in
    db)
        case "${name}" in
            *.syldb|*.syl2db) ;;
            *) echo "ERROR: ${name} is not a sylph database (.syldb or .syl2db)" >&2; exit 1 ;;
        esac
        ;;
    taxonomy)
        case "${name}" in
            *.tsv.gz) ;;
            *) echo "ERROR: ${name} is not a gzipped sylph-tax metadata file (.tsv.gz)" >&2; exit 1 ;;
        esac
        ;;
    *)
        echo "${usage}" >&2
        exit 1
        ;;
esac

mkdir -p sylph_db
wget --tries=5 --continue --waitretry=10 --no-verbose -O "sylph_db/${name}" "${url}"

if [ ! -s "sylph_db/${name}" ]; then
    echo "ERROR: sylph_db/${name} is empty after download" >&2
    exit 1
fi

if [ "${kind}" = "db" ]; then
    # sylph databases are binary (NUL bytes within the first 4 KB of both the
    # .syldb and .syl2db formats); a server error or maintenance page saved
    # under the database name is text and would otherwise be published
    if [ "$(head -c 4096 "sylph_db/${name}" | tr -dc '\000' | wc -c)" -eq 0 ]; then
        echo "ERROR: sylph_db/${name} is text, not a sylph database - the server returned a page instead of the file" >&2
        exit 1
    fi
    # the two-stage format starts with the magic bytes SY2D
    if [ "${name##*.}" = "syl2db" ] && [ "$(head -c 4 "sylph_db/${name}")" != "SY2D" ]; then
        echo "ERROR: sylph_db/${name} does not start with the SY2D magic of a two-stage sylph database" >&2
        exit 1
    fi
fi

if [ "${kind}" = "taxonomy" ]; then
    if ! gzip -t "sylph_db/${name}" 2>/dev/null; then
        echo "ERROR: sylph_db/${name} is not a valid gzip file - corrupt or truncated download" >&2
        exit 1
    fi
    # sylph-tax metadata: genome accession <TAB> d__...;...;s__... lineage
    # (|| true: gzip gets SIGPIPE when head stops reading, fatal under pipefail)
    first="$(gzip -dc "sylph_db/${name}" | head -n 1 || true)"
    if ! printf '%s\n' "${first}" | awk -F'\t' '
            NF >= 2 && $2 ~ /__/ { ok = 1 } END { exit !ok }'; then
        echo "ERROR: sylph_db/${name} is not a sylph-tax metadata file (accession <TAB> d__...;s__... lineage)" >&2
        exit 1
    fi
fi
