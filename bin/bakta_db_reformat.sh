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
#
# The bundled AMRFinderPlus DB is verified present (missing = corrupt or
# truncated download) and then REFRESHED with amrfinder_update: the official
# v6.0 tarball bundles a 2024-era AMRFinderPlus DB that the AMRFinderPlus
# >= 4.x shipped in the pinned bakta biocontainer refuses at runtime
# ("Software requires database version at least 2025-09-22.2"), so without
# the refresh every annotation dies at the AMR expert step (found during the
# T7 acceptance run, 2026-10-06). The refresh (~300 MB from NCBI) runs inside
# the bakta container, which is why FORMAT_BAKTA_DB uses that container.
#
# A DB_VERSION file records the flavor, source URL and the refreshed
# AMRFinderPlus DB version on ONE line for provenance (BAKTA_BAKTA cats it
# into versions.yml, so it must stay single-line; absent in a custom DB, the
# module records 'custom') - the spec requires the light-DB choice to be
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
# corrupt or truncated download (checked before the refresh so a broken
# download is never "self-healed" into a half-valid database)
if [ -z "$(ls -A bakta_db/amrfinderplus-db 2>/dev/null)" ]; then
    echo "ERROR: bakta_db/amrfinderplus-db missing or empty after extraction (corrupt or truncated download?)" >&2
    exit 1
fi

# Refresh the AMRFinderPlus DB to one the container's AMRFinderPlus accepts
# (downloads a new dated subdir and repoints the 'latest' symlink; the
# bundled dir is left in place)
amrfinder_update --force_update --database bakta_db/amrfinderplus-db

amr_version="$(basename "$(readlink -f bakta_db/amrfinderplus-db/latest 2>/dev/null)")"
if [ -z "${amr_version}" ] || [ ! -s "bakta_db/amrfinderplus-db/latest/version.txt" ]; then
    echo "ERROR: amrfinder_update did not produce a valid bakta_db/amrfinderplus-db/latest database" >&2
    exit 1
fi

# Record flavor (full/light, from version.json "type"), source URL and the
# refreshed AMRFinderPlus DB version - single line (see header)
dbtype="$(sed -n 's/.*"type": *"\([^"]*\)".*/\1/p' bakta_db/version.json | head -1)"
printf '%s\n' "db-${dbtype:-unknown} (${url}; amrfinderplus-db ${amr_version})" > bakta_db/DB_VERSION

if [ ! -s bakta_db/DB_VERSION ]; then
    echo "ERROR: expected file 'bakta_db/DB_VERSION' missing or empty" >&2
    exit 1
fi
