#!/bin/bash
#
# rgi_load_kmers.sh <card_json> <wildcard_dir> <kmer_size>
#
# Loads the WildCARD k-mer database into the RGI --local database (localDB/
# in the current directory), for rgi kmer_query (RGI_KMER, pathogen of
# origin). kmer_query reads localDB/<k>mer_database.json and
# localDB/amr_<k>mer.txt; when they are missing it logs an error and exits 0,
# so without this step every pathogen-of-origin result was empty.
# The WildCARD variants download ships them as <k>_kmer_db.json and
# all_amr_<k>mers.txt (gzipped in the tarball). rgi load requires --card_json
# on every call; it adds to the localDB/loaded_databases.json of the earlier
# CARD/WildCARD load instead of replacing it.
#
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: rgi_load_kmers.sh <card_json> <wildcard_dir> <kmer_size>" >&2
    exit 2
fi
card_json="$1"
wildcard_dir="$2"
k="$3"

kmer_db="${wildcard_dir}/${k}_kmer_db.json"
amr_kmers="${wildcard_dir}/all_amr_${k}mers.txt"
for f in "${kmer_db}" "${amr_kmers}"; do
    if [ ! -f "${f}" ] && [ -f "${f}.gz" ]; then
        gunzip "${f}.gz"
    fi
    if [ ! -f "${f}" ]; then
        echo "ERROR: ${f} not found: the WildCARD directory has no ${k}-mer database," >&2
        echo "       which rgi kmer_query (pathogen of origin) needs. Check --rgi_kmer_size" >&2
        echo "       against the *_kmer_db.json files in the WildCARD download." >&2
        exit 1
    fi
done

rgi load --card_json "${card_json}" --kmer_database "${kmer_db}" --amr_kmers "${amr_kmers}" --kmer_size "${k}" --local

for f in "localDB/${k}mer_database.json" "localDB/amr_${k}mer.txt"; do
    if [ ! -s "${f}" ]; then
        echo "ERROR: rgi load did not create ${f}" >&2
        exit 1
    fi
done
echo "Loaded the ${k}-mer database into localDB/"
