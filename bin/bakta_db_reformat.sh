#!/bin/bash
#
# Download and unpack the Bakta database for MAG annotation (design doc
# Section 4.5, T7).
# Usage: bakta_db_reformat.sh <bakta-db-tarball-url>
#   (e.g. https://zenodo.org/record/14916843/files/db.tar.xz       - full, 31.9 GB
#      or https://zenodo.org/record/14916843/files/db-light.tar.xz - light, 1.3 GB)
#
# Produces bakta_db/ with the schema-6 layout Bakta 1.12.x expects for --db
# (version.json, amrfinderplus-db/, the pre-annotated sequence databases).
# --strip-components=1 normalizes the tarball's top directory (db/ for full,
# db-light/ for light) so both flavors land in the same bakta_db/ path.
# The bundled AMRFinderPlus DB is verified, never rebuilt: a download whose
# amrfinderplus-db is missing is corrupt and must fail, not self-heal via
# amrfinder_update. A DB_VERSION file records the flavor and source URL for
# provenance (BAKTA_BAKTA reads it into versions.yml; absent in a custom DB,
# the module records 'custom') - the spec requires the light-DB choice to be
# recorded because it changes results.
#
set -euo pipefail

url="${1:?usage: bakta_db_reformat.sh <bakta-db-tarball-url>}"
tarball="${url##*/}"

wget --tries=5 --continue --waitretry=10 --no-verbose -O "${tarball}" "${url}"

mkdir -p bakta_db
tar -xJf "${tarball}" --strip-components=1 -C bakta_db
rm -f "${tarball}"

if [ ! -s bakta_db/version.json ]; then
    echo "ERROR: bakta_db/version.json missing or empty after extraction" >&2
    exit 1
fi

# Bakta 1.12.x requires database schema major version 6
major="$(sed -n 's/.*"major": *\([0-9][0-9]*\).*/\1/p' bakta_db/version.json | head -1)"
if [ "${major}" != "6" ]; then
    echo "ERROR: bakta_db/version.json reports schema major '${major:-missing}', expected 6 (required by Bakta 1.12.x)" >&2
    exit 1
fi

# The tarball bundles the AMRFinderPlus DB; a missing/empty one means a
# corrupt or truncated download
if [ -z "$(ls -A bakta_db/amrfinderplus-db 2>/dev/null)" ]; then
    echo "ERROR: bakta_db/amrfinderplus-db missing or empty after extraction (corrupt or truncated download?)" >&2
    exit 1
fi

# Record flavor (full/light, from version.json "type") plus source URL
dbtype="$(sed -n 's/.*"type": *"\([^"]*\)".*/\1/p' bakta_db/version.json | head -1)"
printf '%s\n' "db-${dbtype:-unknown} (${url})" > bakta_db/DB_VERSION

if [ ! -s bakta_db/DB_VERSION ]; then
    echo "ERROR: expected file 'bakta_db/DB_VERSION' missing or empty" >&2
    exit 1
fi
