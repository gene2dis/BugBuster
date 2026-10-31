#!/bin/bash
#
# Build the CARD reference K-MARVEL reads, as upstream's get_ref_db.sh does.
# Usage: kmarvel_db_reformat.sh <card-data-tarball-url> <md5> <dbversion>
#        (e.g. https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2)
#
# Produces kmarvel_db/ with
#   CARD_protein_all.fasta   protein homolog + overexpression + variant model
#                            FASTAs, concatenated (the --reference)
#   metadata.tsv             CARD aro_index.tsv (the --metadata)
#   snps.txt                 CARD snps.txt (the --snpdata; the copy in
#                            K-MARVEL's repo db/ is this same file)
#   DB_VERSION               provenance, read into the KMARVEL versions.yml
#   CARD_protein_all.fasta.p*  its BLAST index (--reference is also the
#                            blastp -db), built once here; KMARVEL links it.
#                            A --custom_kmarvel_db without one gets the index
#                            built per task instead
# (~2.5 MB FASTA).
#
set -euo pipefail

card_url="${1:?usage: kmarvel_db_reformat.sh <card-data-tarball-url> <md5> <dbversion>}"
card_md5="${2:?usage: kmarvel_db_reformat.sh <card-data-tarball-url> <md5> <dbversion>}"
dbversion="${3:?usage: kmarvel_db_reformat.sh <card-data-tarball-url> <md5> <dbversion>}"

wget --tries=5 --continue --waitretry=10 --no-verbose -O card.tar.bz2 "${card_url}"
echo "${card_md5}  card.tar.bz2" | md5sum -c -

mkdir -p card_data kmarvel_db
tar -xjf card.tar.bz2 -C card_data

for f in protein_fasta_protein_homolog_model.fasta \
         protein_fasta_protein_overexpression_model.fasta \
         protein_fasta_protein_variant_model.fasta \
         aro_index.tsv snps.txt; do
    if [ ! -s "card_data/${f}" ]; then
        echo "ERROR: '${f}' missing or empty in the CARD data tarball" >&2
        exit 1
    fi
done

cat card_data/protein_fasta_protein_homolog_model.fasta \
    card_data/protein_fasta_protein_overexpression_model.fasta \
    card_data/protein_fasta_protein_variant_model.fasta \
    > kmarvel_db/CARD_protein_all.fasta
mv card_data/aro_index.tsv kmarvel_db/metadata.tsv
mv card_data/snps.txt kmarvel_db/snps.txt
printf '%s\n' "${dbversion}" > kmarvel_db/DB_VERSION

rm -rf card_data card.tar.bz2

if ! grep -q '^>' kmarvel_db/CARD_protein_all.fasta; then
    echo "ERROR: kmarvel_db/CARD_protein_all.fasta has no FASTA records" >&2
    exit 1
fi

makeblastdb -in kmarvel_db/CARD_protein_all.fasta -dbtype prot -blastdb_version 5 \
    -title CARD_protein_all -out kmarvel_db/CARD_protein_all.fasta > /dev/null
