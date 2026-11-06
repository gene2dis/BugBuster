#!/bin/bash
#
# Download, verify and unpack a GUNC DIAMOND reference (proGenomes 2.1: a
# 7.2 GB download, 13.5 GB unpacked), the file GUNC 1.1.1's download_db
# fetches, checked against both md5s it records. GUNC reads which reference
# it is from the .dmnd file name, so the file keeps it.
# Usage: gunc_db_reformat.sh <dmnd.gz-url> <gz-md5> <dmnd-file-name> <dmnd-md5>
#
set -euo pipefail

usage="usage: gunc_db_reformat.sh <dmnd.gz-url> <gz-md5> <dmnd-file-name> <dmnd-md5>"
db_url="${1:?${usage}}"
gz_md5="${2:?${usage}}"
db_file="${3:?${usage}}"
db_md5="${4:?${usage}}"

wget --tries=5 --continue --waitretry=10 --no-verbose -O "${db_file}.gz" "${db_url}"
echo "${gz_md5}  ${db_file}.gz" | md5sum -c -
gunzip "${db_file}.gz"
echo "${db_md5}  ${db_file}" | md5sum -c -
