#!/bin/bash
#
# Print the CARD version K-MARVEL ran against, for the "card:" entry of the
# KMARVEL versions.yml.
# Usage: kmarvel_db_version.sh <kmarvel_db_dir>
#
#   - DB_VERSION present and not blank (FORMAT_KMARVEL_DB builds, or a copy
#     of one): its first line
#   - otherwise (a --custom_kmarvel_db without DB_VERSION):
#     "custom (CARD_protein_all.fasta md5 <md5>)", so the reference can
#     still be traced
#
# Double quotes, backslashes and carriage returns are dropped: the value is
# written inside a double-quoted YAML string.
#
set -euo pipefail

db_dir="${1:?usage: kmarvel_db_version.sh <kmarvel_db_dir>}"

version=""
if [ -s "${db_dir}/DB_VERSION" ]; then
    version=$(head -n 1 "${db_dir}/DB_VERSION" | tr -d '"\\\r')
fi

if [ -z "${version//[[:space:]]/}" ]; then
    md5=$(md5sum "${db_dir}/CARD_protein_all.fasta" | cut -d' ' -f1)
    version="custom (CARD_protein_all.fasta md5 ${md5})"
fi

printf '%s\n' "${version}"
