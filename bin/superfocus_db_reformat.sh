#!/bin/bash
#
# Provision the SUPER-FOCUS DB_90 database for the read-level functional
# branch (design doc Section 4.6.3, task T8b, Q17).
# Usage: superfocus_db_reformat.sh <aligner> <zip-url> <zip-md5> <pks-url> <pks-md5> <provenance>
#   aligner     diamond | mmseqs2 (--superfocus_aligner)
#   zip-url     figshare download of the prebuilt DB_90 archive for that
#               aligner (diamond format v3 90_clusters.db.dmnd.zip, or
#               mmseqs_90.zip) - CC0, open.flinders.edu.au
#   zip-md5     the md5 figshare publishes for that archive
#   pks-url     database_PKs.txt pinned to the SUPER-FOCUS v1.8 tag (the 1.8
#               package ships no db/ folder, so the subsystem table is
#               fetched separately)
#   pks-md5     its md5
#   provenance  human-readable release string for DB_VERSION
#
# Produces superfocus_db/, the database ROOT `superfocus -b` expects:
#   superfocus_db/db/database_PKs.txt
#   superfocus_db/db/static/diamond/90_clusters.db.dmnd     (diamond), or
#   superfocus_db/db/static/mmseqs2/90_clusters.db*         (mmseqs2)
#   superfocus_db/DB_VERSION   single line (read by SUPERFOCUS into versions.yml)
# Only the selected aligner's files are downloaded (the formats are
# aligner-specific). Archive and table are md5-verified before use; the
# unpacked database must contain the aligner's main file, non-empty.
#
set -euo pipefail

usage="usage: superfocus_db_reformat.sh <diamond|mmseqs2> <zip-url> <zip-md5> <pks-url> <pks-md5> <provenance>"
aligner="${1:?${usage}}"
zip_url="${2:?${usage}}"
zip_md5="${3:?${usage}}"
pks_url="${4:?${usage}}"
pks_md5="${5:?${usage}}"
provenance="${6:?${usage}}"

case "${aligner}" in
    diamond) static_dir=diamond; main_file=90_clusters.db.dmnd ;;
    mmseqs2) static_dir=mmseqs2; main_file=90_clusters.db ;;
    *) echo "ERROR: unknown aligner '${aligner}' (expected diamond or mmseqs2)" >&2; exit 1 ;;
esac

fetch() {  # fetch <url> <output>
    wget --tries=5 --continue --waitretry=10 --no-verbose -O "$2" "$1"
}

verify_md5() {  # verify_md5 <file> <expected> <label>
    local actual
    actual="$(md5sum "$1" | cut -d' ' -f1)"
    if [ "${actual}" != "$2" ]; then
        echo "ERROR: md5 mismatch for $3 (got ${actual}, expected $2) - corrupt or truncated download" >&2
        exit 1
    fi
}

target="superfocus_db/db/static/${static_dir}"
mkdir -p "${target}"

fetch "${pks_url}" superfocus_db/db/database_PKs.txt
verify_md5 superfocus_db/db/database_PKs.txt "${pks_md5}" database_PKs.txt
if [ "$(head -1 superfocus_db/db/database_PKs.txt)" != "$(printf 'Primary key\tLevel 1\tLevel 2\tLevel 3/Subsystem Name')" ]; then
    echo "ERROR: database_PKs.txt does not have the SUPER-FOCUS subsystem table header" >&2
    exit 1
fi

archive="${zip_url##*/}.zip"
fetch "${zip_url}" "${archive}"
verify_md5 "${archive}" "${zip_md5}" "the ${aligner} DB_90 archive"
# -j: the archives' internal folder layout is not relied on; every member
# lands directly in db/static/<aligner>/
unzip -q -o -j "${archive}" -d "${target}"
rm -f "${archive}"

if [ ! -s "${target}/${main_file}" ]; then
    echo "ERROR: ${target}/${main_file} missing or empty after unpacking - not a SUPER-FOCUS ${aligner} DB_90 archive" >&2
    exit 1
fi
if [ "${aligner}" = "mmseqs2" ] && [ ! -s "${target}/${main_file}.dbtype" ]; then
    echo "ERROR: ${target}/${main_file}.dbtype missing - incomplete MMseqs2 database" >&2
    exit 1
fi

printf '%s\n' "${provenance}; ${aligner} DB_90 (${zip_url}, md5 ${zip_md5})" > superfocus_db/DB_VERSION
