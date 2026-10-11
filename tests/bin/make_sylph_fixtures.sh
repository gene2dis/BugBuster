#!/bin/bash
#
# Regenerate the sylph module fixtures in tests/data/sylph/ (used by
# tests/modules/sylph_profile.nf.test and sylph_taxprof.nf.test):
#   sylph_test_db.syldb       sylph 1.0.0 database of two genomes from
#                             tests/data/decontamination: phage lambda
#                             (host_mock.fasta) as GCF_000840245.1 and phiX174
#                             as GCF_000819615.1 (default -c 200)
#   sylph_test_R{1,2}.fastq.gz  1,000 error-light read pairs simulated from
#                             lambda (~6x) + 500 random pairs (unknown reads)
#   sylph_tax_fixture.tsv.gz  sylph-tax metadata for the two accessions
#   sylph_profile_hit.tsv     real sylph profile -u of those reads (lambda)
#   sylph_profile_empty.tsv   real header-only profile of the random pairs
#
# tests/data/taxonomy/sampleA_sylph_gtdb.sylphmpa is NOT made here: it is a
# verbatim sylph-tax 1.9.2 capture of a real sylph 1.0.0 profile against the
# 20 GB GTDB r232 .syl2db (simulated E. coli K-12 + B. subtilis 168 reads),
# kept as committed data.
#
# Reads and metadata are byte-stable (fixed seed, gzip mtime 0). Requires
# docker and python3; no network.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
OUT_DIR="${REPO_DIR}/tests/data/sylph"
SYLPH_IMG=$(grep -o "'quay.io/biocontainers/sylph:[^']*'" "${REPO_DIR}/modules/local/sylph_profile/main.nf" | tr -d "'")
mkdir -p "${OUT_DIR}"

WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT

cp "${REPO_DIR}/tests/data/decontamination/host_mock.fasta" "${WORK}/GCF_000840245.1_genomic.fna"
cp "${REPO_DIR}/tests/data/decontamination/phiX174.fasta"   "${WORK}/GCF_000819615.1_genomic.fna"

python3 -I - "${WORK}/GCF_000840245.1_genomic.fna" "${WORK}" <<'EOF'
import gzip, random, sys
genome_fa, out = sys.argv[1], sys.argv[2]
random.seed(42)
seq = ''.join(l.strip().upper() for l in open(genome_fa) if not l.startswith('>'))
comp = str.maketrans('ACGTN', 'TGCAN')
rc = lambda s: s.translate(comp)[::-1]
mut = lambda s: ''.join(random.choice('ACGT') if random.random() < 0.001 else c for c in s)
L, INS = 150, 300
recs = []
for i in range(1000):
    p = random.randint(0, len(seq) - INS)
    frag = seq[p:p + INS]
    if random.random() < 0.5:
        frag = rc(frag)
    recs.append((f'lambda_{i}', mut(frag[:L]), mut(rc(frag)[:L])))
rnd = lambda: ''.join(random.choice('ACGT') for _ in range(L))
unknown = [(f'unknown_{i}', rnd(), rnd()) for i in range(500)]
random.shuffle(recs := recs + unknown)
q = 'I' * L
def write(path, records, mate):
    with open(path, 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as gz:
        for name, r1, r2 in records:
            gz.write(f'@{name}/{mate}\n{r1 if mate == 1 else r2}\n+\n{q}\n'.encode())
for mate in (1, 2):
    write(f'{out}/sylph_test_R{mate}.fastq.gz', recs, mate)
    write(f'{out}/unknown_R{mate}.fastq.gz', unknown, mate)
EOF

printf 'GCF_000840245.1\td__Viruses;p__Uroviricota;c__Caudoviricetes;o__Caudovirales;f__Lambdaviridae;g__Lambdavirus;s__Lambdavirus lambda\nGCF_000819615.1\td__Viruses;p__Phixviricota;c__Malgrandaviricetes;o__Petitvirales;f__Microviridae;g__Sinsheimervirus;s__Sinsheimervirus phiX174\n' \
    | gzip -n > "${WORK}/sylph_tax_fixture.tsv.gz"

docker run --rm -u "$(id -u):$(id -g)" -v "${WORK}:/w" -w /w "${SYLPH_IMG}" bash -c '
    set -euo pipefail
    sylph sketch GCF_000840245.1_genomic.fna GCF_000819615.1_genomic.fna -o sylph_test_db 2>/dev/null
    sylph profile sylph_test_db.syldb -1 sylph_test_R1.fastq.gz -2 sylph_test_R2.fastq.gz -u -t 1 -o sylph_profile_hit.tsv 2>/dev/null
    sylph profile sylph_test_db.syldb -1 unknown_R1.fastq.gz -2 unknown_R2.fastq.gz -u -t 1 -o sylph_profile_empty.tsv 2>/dev/null
'

if ! grep -q GCF_000840245.1 "${WORK}/sylph_profile_hit.tsv"; then
    echo "ERROR: lambda not detected in sylph_profile_hit.tsv" >&2
    exit 1
fi
if [ "$(wc -l < "${WORK}/sylph_profile_empty.tsv")" -ne 1 ]; then
    echo "ERROR: sylph_profile_empty.tsv is not header-only" >&2
    exit 1
fi

for f in sylph_test_db.syldb sylph_test_R1.fastq.gz sylph_test_R2.fastq.gz \
         sylph_tax_fixture.tsv.gz sylph_profile_hit.tsv sylph_profile_empty.tsv; do
    cp "${WORK}/${f}" "${OUT_DIR}/${f}"
    echo "wrote ${OUT_DIR}/${f}"
done
