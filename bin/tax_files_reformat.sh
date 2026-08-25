#!/bin/bash
#
# Download the NCBI taxdump and keep nodes.dmp/names.dmp.
# Usage: tax_files_reformat.sh <taxdump.tar.gz-url>
#
set -euo pipefail

wget --tries=5 --continue --waitretry=10 -O taxdump.tar.gz \
    "${1:?usage: tax_files_reformat.sh <taxdump.tar.gz-url>}"

mkdir tax_files
tar -xzf taxdump.tar.gz

if [ ! -f nodes.dmp ] || [ ! -f names.dmp ]; then
    echo "ERROR: nodes.dmp/names.dmp missing after extracting taxdump.tar.gz" >&2
    exit 1
fi
mv nodes.dmp names.dmp tax_files
rm -f ./*.dmp gc.prt readme.txt taxdump.tar.gz
