#!/bin/bash
#
# Download the NCBI COG definitions table used to map eggNOG COG ids to COG
# functional-category letters in the contig functional branch (design doc
# Q16, owner option b).
# Usage: cog_db_reformat.sh <cog-def-table-url>
#   (e.g. https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.def.tab)
#
# Produces <name>.def.tab (e.g. cog-24.def.tab; ~410 KB) in the working
# directory. The release's checksums.md5 (same directory) is fetched and the
# table verified against it; the layout is then checked (tab-separated, COG
# id in column 1, category letters in column 2) so a wrong or truncated file
# fails here instead of in aggregation. The filename carries the release and
# is recorded as provenance by AGGREGATE_FUNCTIONS.
#
set -euo pipefail

url="${1:?usage: cog_db_reformat.sh <cog-def-table-url>}"
name="${url##*/}"
base="${url%/*}"

wget --tries=5 --continue --waitretry=10 --no-verbose -O "${name}" "${url}"
wget --tries=5 --waitretry=10 --no-verbose -O checksums.md5 "${base}/checksums.md5"

expected="$(awk -v f="${name}" '$2 == f || $2 == "*" f { print $1 }' checksums.md5 | head -1)"
actual="$(md5sum "${name}" | cut -d' ' -f1)"
if [ -z "${expected}" ] || [ "${expected}" != "${actual}" ]; then
    echo "ERROR: md5 mismatch for ${name} (got ${actual}, published ${expected:-missing}) - corrupt or truncated download" >&2
    exit 1
fi
rm -f checksums.md5

if ! awk -F'\t' '
        { gsub(/[ \r]+$/, "", $2) }
        $1 !~ /^COG[0-9]+$/ || $2 !~ /^[A-Z]+$/ { bad = 1; exit }
        END { exit (bad || NR == 0) }' "${name}"; then
    echo "ERROR: ${name} is not an NCBI COG definitions table (COG id <TAB> category letters <TAB> ...)" >&2
    exit 1
fi
