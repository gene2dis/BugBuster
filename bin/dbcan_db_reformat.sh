#!/bin/bash
#
# Download the dbCAN database files needed by run_dbcan v5 protein-mode
# CAZyme annotation (the tool's DATABASES_CAZYME set).
# Usage: dbcan_db_reformat.sh <dbcan-release-base-url> ["<local file>:<bytes> ..."]
#   (e.g. https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026; the
#   byte counts come from the registry and are checked after the download)
#
# Produces dbcan_db/ with the layout run_dbcan expects for --db_dir:
#   CAZy.dmnd, dbCAN.hmm, dbCAN-sub.hmm, fam-substrate-mapping.tsv
# (~7.4 GB total). The release serves the sub-HMM file as dbCAN_sub.hmm
# (underscore) but the tool looks for dbCAN-sub.hmm (hyphen), so it is
# renamed on download. No hmmpress step: run_dbcan v5 uses pyhmmer, which
# reads plain .hmm files. A DB_VERSION file records the release string for
# provenance (RUN_DBCAN reads it into versions.yml; absent in a custom DB,
# the module records 'custom').
#
set -euo pipefail

base_url="${1:?usage: dbcan_db_reformat.sh <dbcan-release-base-url> [\"<local file>:<bytes> ...\"]}"
BYTES="${2:-}"

# remote-name:local-name pairs (only the sub-HMM file differs)
files="CAZy.dmnd:CAZy.dmnd dbCAN.hmm:dbCAN.hmm dbCAN_sub.hmm:dbCAN-sub.hmm fam-substrate-mapping.tsv:fam-substrate-mapping.tsv"

mkdir -p dbcan_db
for pair in ${files}; do
    remote="${pair%%:*}"
    local_name="${pair##*:}"
    wget --tries=5 --continue --waitretry=10 --no-verbose -O "dbcan_db/${local_name}" "${base_url}/${remote}"
done

# Byte counts from the registry (optional last argument, "name:bytes ..."):
# the files are served uncompressed, so a download cut at a block boundary or
# replaced by a proxy page has no format to fail on (review R6-13)
for pair in ${BYTES}; do
    f="${pair%:*}"
    want="${pair##*:}"
    got=$(stat -c %s "dbcan_db/${f}" 2>/dev/null || echo missing)
    if [ "${got}" != "${want}" ]; then
        echo "ERROR: dbcan_db/${f} is ${got} bytes, the release has ${want} - a truncated download, or the server now serves another file under this URL (the registry entry pins this release: pass a copy you trust with --custom_dbcan_db, or update config/databases.config with a new key)" >&2
        exit 1
    fi
done

printf '%s\n' "${base_url##*/}" > dbcan_db/DB_VERSION

# Verify the database landed before finishing, so exit status reflects the DB
for f in CAZy.dmnd dbCAN.hmm dbCAN-sub.hmm fam-substrate-mapping.tsv DB_VERSION; do
    if [ ! -s "dbcan_db/${f}" ]; then
        echo "ERROR: expected file 'dbcan_db/${f}' missing or empty after download" >&2
        exit 1
    fi
done
