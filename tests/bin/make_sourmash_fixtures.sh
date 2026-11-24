#!/bin/bash
#
# Regenerate the sourmash module fixtures in tests/data/sourmash/ (used by
# tests/modules/sourmash.nf.test; roadmap Phase 25, RV25):
#   two_genome_k31.zip        sourmash 4.9.4 database (k=31, scaled=1000, abund)
#                             of two made-up genomes: GCF_900000001.1 (120 kb)
#                             and GCF_900000002.1 (100 kb), seeded random DNA.
#                             Both are above sourmash gather's default 50 kbp
#                             threshold (phage lambda, 48.5 kb, is not)
#   two_genome_lineages.csv   their lineages, in the GTDB lineages CSV layout
#                             SOURMASH_TAX_PREPARE reads
#   sm_test_R{1,2}.fastq.gz   1,000 error-free 150 bp pairs simulated from
#                             GCF_900000001.1 (~2.5x)
# The no-match case uses tests/data/reads/test_sample1 (B. fragilis reads).
#
# stub_db.zip and stub_lineages.csv in the same folder are NOT made here: they
# are hand-written placeholders for the -stub tests (sourmash is not run).
#
# Reads and lineages are byte-stable (fixed seed, gzip mtime 0); the zip holds
# sourmash's own timestamps, so a rebuild differs in bytes but not in content.
# Requires docker and python3; no network.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
OUT_DIR="${REPO_DIR}/tests/data/sourmash"
SOURMASH_IMG=$(grep -o "'quay.io/biocontainers/sourmash:[^']*'" "${REPO_DIR}/modules/local/sourmash/main.nf" | tr -d "'")
if [ -z "${SOURMASH_IMG}" ]; then
    echo "ERROR: could not read the sourmash image from modules/local/sourmash/main.nf" >&2
    exit 1
fi
mkdir -p "${OUT_DIR}"

WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT

python3 -I - "${WORK}" "${OUT_DIR}" <<'EOF'
import gzip, random, sys
work, out = sys.argv[1], sys.argv[2]
rng = random.Random(25)
genomes = {"GCF_900000001.1": 120000, "GCF_900000002.1": 100000}
seqs = {}
for acc, n in genomes.items():
    seqs[acc] = "".join(rng.choice("ACGT") for _ in range(n))
    with open(f"{work}/{acc}.fna", "w") as f:
        f.write(f">{acc} synthetic genome\n")
        f.writelines(seqs[acc][i:i + 80] + "\n" for i in range(0, n, 80))
with open(f"{out}/two_genome_lineages.csv", "w") as f:
    f.write("ident,superkingdom,phylum,class,order,family,genus,species\n")
    f.write("GCF_900000001.1,d__Bacteria,p__Synthetica,c__Synthetia,o__Synthetales,f__Synthetaceae,g__Synthetus,s__Synthetus alpha\n")
    f.write("GCF_900000002.1,d__Bacteria,p__Synthetica,c__Synthetia,o__Synthetales,f__Synthetaceae,g__Synthetus,s__Synthetus beta\n")
comp = str.maketrans("ACGT", "TGCA")
g = seqs["GCF_900000001.1"]
r1, r2 = [], []
for k in range(1000):
    st = rng.randrange(0, len(g) - 300)
    frag = g[st:st + 300]
    if rng.random() < 0.5:
        frag = frag.translate(comp)[::-1]
    r1.append(f"@sm{k}/1\n{frag[:150]}\n+\n{'I' * 150}\n")
    r2.append(f"@sm{k}/2\n{frag.translate(comp)[::-1][:150]}\n+\n{'I' * 150}\n")
for name, recs in (("sm_test_R1.fastq.gz", r1), ("sm_test_R2.fastq.gz", r2)):
    with open(f"{out}/{name}", "wb") as fh:
        with gzip.GzipFile(fileobj=fh, mode="wb", mtime=0) as gz:
            gz.write("".join(recs).encode())
EOF

docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${WORK}:/w" -w /w "${SOURMASH_IMG}" bash -euo pipefail -c '
    for acc in GCF_900000001.1 GCF_900000002.1; do
        sourmash sketch dna -p k=31,scaled=1000,abund "${acc}.fna" --name "${acc}" -o "${acc}.sig" 2> /dev/null
    done
    sourmash sig cat GCF_900000001.1.sig GCF_900000002.1.sig -o two_genome_k31.zip 2> /dev/null
'
cp "${WORK}/two_genome_k31.zip" "${OUT_DIR}/"
ls -l "${OUT_DIR}"
