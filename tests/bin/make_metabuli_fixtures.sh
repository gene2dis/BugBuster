#!/bin/bash
#
# Regenerate the Metabuli module fixtures in tests/data/metabuli/ (used by
# tests/modules/metabuli_classify.nf.test):
#   mini_db/            Metabuli 1.2.0 database of two genomes from
#                       tests/data/decontamination: phage lambda
#                       (host_mock.fasta, taxid 10710) and phiX174 (10847),
#                       with its taxdump in mini_db/ncbi-taxdump/ (2025 NCBI
#                       style: Viruses is an 'acellular root', no 'domain')
#   contigs.fa          three MEGAHIT-style contigs: 3 kb of lambda, 2.5 kb
#                       of phiX174 and 1.5 kb of seeded random sequence (no
#                       call)
#
# taxonomyDB is not byte-stable between builds (Metabuli serialises it in a
# varying order); the calls on contigs.fa are. Requires docker and python3;
# no network.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
OUT_DIR="${REPO_DIR}/tests/data/metabuli"
METABULI_IMG=$(grep -o "'quay.io/biocontainers/metabuli:[^']*'" "${REPO_DIR}/modules/local/metabuli_classify/main.nf" | tr -d "'")
DECONT="${REPO_DIR}/tests/data/decontamination"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

mkdir -p "${TMP_DIR}/tax"
nodes_row() { printf '%s\t|\t%s\t|\t%s\t|\t\t|\t9\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\n' "$1" "$2" "$3"; }
names_row() { printf '%s\t|\t%s\t|\t\t|\tscientific name\t|\n' "$1" "$2"; }
{
    nodes_row 1 1 'no rank'
    nodes_row 10239 1 'acellular root'
    nodes_row 2731619 10239 class
    nodes_row 10710 2731619 species
    nodes_row 2843429 10239 class
    nodes_row 10847 2843429 species
} > "${TMP_DIR}/tax/nodes.dmp"
{
    names_row 1 root
    names_row 10239 Viruses
    names_row 2731619 Caudoviricetes
    names_row 10710 'Lambdavirus lambda'
    names_row 2843429 Malgrandaviricetes
    names_row 10847 'Sinsheimervirus phiX174'
} > "${TMP_DIR}/tax/names.dmp"
: > "${TMP_DIR}/tax/merged.dmp"
: > "${TMP_DIR}/tax/delnodes.dmp"
printf 'accession\taccession.version\ttaxid\tgi\nNC_001416\tNC_001416.1\t10710\t0\nNC_001422\tNC_001422.1\t10847\t0\n' \
    > "${TMP_DIR}/acc2taxid.map"
cat "${DECONT}/host_mock.fasta" "${DECONT}/phiX174.fasta" > "${TMP_DIR}/genomes.fa"
echo "${TMP_DIR}/genomes.fa" > "${TMP_DIR}/fasta_list.txt"

docker run --rm -u "$(id -u):$(id -g)" -v "${TMP_DIR}:${TMP_DIR}" -w "${TMP_DIR}" "${METABULI_IMG}" \
    metabuli build mini_db fasta_list.txt acc2taxid.map --taxonomy-path tax --threads 2 --max-ram 2 \
    --split-num 8 --db-date 2026-10-29 --db-name bugbuster_test > "${TMP_DIR}/build.log" 2>&1

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"
mv "${TMP_DIR}/mini_db" "${OUT_DIR}/mini_db"
mkdir -p "${OUT_DIR}/mini_db/ncbi-taxdump"
cp "${TMP_DIR}/tax/nodes.dmp" "${TMP_DIR}/tax/names.dmp" "${OUT_DIR}/mini_db/ncbi-taxdump/"

python3 - "${DECONT}" "${OUT_DIR}/contigs.fa" <<'PYEOF'
import random
import sys

decont, out = sys.argv[1], sys.argv[2]


def seq(path):
    return ''.join(line.strip() for line in open(path) if not line.startswith('>'))


random.seed(7)
lam = seq(f'{decont}/host_mock.fasta')
phx = seq(f'{decont}/phiX174.fasta')
rnd = ''.join(random.choice('ACGT') for _ in range(1500))
with open(out, 'w') as fh:
    fh.write('>k141_1 flag=1 multi=5.0 len=3000\n' + lam[10000:13000] + '\n')
    fh.write('>k141_2 flag=1 multi=4.0 len=2500\n' + phx[1000:3500] + '\n')
    fh.write('>k141_3 flag=0 multi=2.0 len=1500\n' + rnd + '\n')
PYEOF

echo "Wrote ${OUT_DIR}:"
ls -la "${OUT_DIR}" "${OUT_DIR}/mini_db"
