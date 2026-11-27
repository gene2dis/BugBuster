#!/bin/bash
#
# Regenerates tests/data/semibin/ (roadmap Phase 29, CL3): contig sets on which
# SemiBin2 2.6.0 finds no bins, each with a BAM of reads tiled along its
# contigs, for tests/modules/semibin.nf.test.
#   short_contigs.fa / .bam  two contigs of 1,200 and 1,500 bp: SemiBin2 stops
#                            with "all are shorter than 2500 basepairs"
#   mid_contigs.fa / .bam    four contigs of 2,800-3,400 bp: SemiBin2 stops
#                            with "no must-link pairs can be generated"
# Sequences are random (seeded), not real data. Reads are 100 bp, every 5 bp
# on the forward strand. The BAMs are sorted with the samtools of the SEMIBIN
# image (--no-PG: no command line in the header), so a rerun gives the same
# bytes. Needs python3 and docker.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUT="$( cd "${SCRIPT_DIR}/../data" && pwd )/semibin"
IMAGE="quay.io/biocontainers/semibin:2.6.0--pyhdfd78af_0"
mkdir -p "${OUT}"

python3 - "${OUT}" <<'PYEOF'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
sets = {"short_contigs": [1200, 1500], "mid_contigs": [3000, 3200, 3400, 2800]}
rng = random.Random(29)
for name, lengths in sets.items():
    seqs = {f"contig_{i}": "".join(rng.choice("ACGT") for _ in range(n)) for i, n in enumerate(lengths, 1)}
    with open(out / f"{name}.fa", "w") as f:
        for cid, s in seqs.items():
            f.write(f">{cid}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    with open(out / f"{name}.sam", "w") as f:
        for cid, s in seqs.items():
            f.write(f"@SQ\tSN:{cid}\tLN:{len(s)}\n")
        k = 0
        for cid, s in seqs.items():
            for p in range(0, len(s) - 100, 5):
                k += 1
                f.write(f"r{k}\t0\t{cid}\t{p + 1}\t60\t100M\t*\t0\t0\t{s[p:p + 100]}\t{'I' * 100}\n")
PYEOF

docker run --rm -u "$(id -u):$(id -g)" -v "${OUT}:/d" -w /d "${IMAGE}" bash -c '
    for n in short_contigs mid_contigs; do
        samtools sort --no-PG -o ${n}.bam ${n}.sam
        rm ${n}.sam
    done'
ls -l "${OUT}"
