#!/usr/bin/env python3
"""Generate the tiny deterministic FASTQ fixtures used by the singleton QC
regression test (tests/main.nf.test, audit items #2/#9/#16).

Reads tests/data/decontamination/host_mock.fasta (lambda phage, the mock host
genome) and writes tests/data/reads/singleton_test_{R1,R2,S}.fastq.gz:

- R1/R2: 6 non-host pairs + 4 host-derived pairs (proper FR orientation,
  ~300 bp insert). The host pairs must be removed by decontamination,
  leaving exactly 6 clean pairs (12 reads).
- S: 4 non-host singletons + 6 host-derived singletons, leaving exactly
  4 clean singletons.

Reads are exact 100-mers at uniform Q40 with no long G runs, so fastp
(default pipeline settings) keeps all of them and the expected clean counts
are deterministic. Non-host sequence comes from a fixed-seed random generator;
a chance 21 bp+ local alignment against lambda/phiX is vanishingly unlikely.
"""

import gzip
import random
from pathlib import Path

READ_LEN = 100
INSERT = 300
QUAL = "I" * READ_LEN  # Q40

BASE = Path(__file__).resolve().parent.parent
HOST_FASTA = BASE / "data" / "decontamination" / "host_mock.fasta"
OUT_DIR = BASE / "data" / "reads"

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def random_seq(rng, length):
    while True:
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        if "G" * 10 not in seq:  # keep clear of fastp --trim_poly_g
            return seq


def fastq_record(name, seq):
    return f"@{name}\n{seq}\n+\n{QUAL}\n"


def main():
    host = "".join(
        line.strip() for line in HOST_FASTA.read_text().splitlines()
        if not line.startswith(">")
    )
    rng = random.Random(42)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    r1, r2, s = [], [], []

    # 6 non-host pairs (survive decontamination)
    for i in range(6):
        frag = random_seq(rng, INSERT)
        r1.append(fastq_record(f"nonhost_pair_{i}/1", frag[:READ_LEN]))
        r2.append(fastq_record(f"nonhost_pair_{i}/2", revcomp(frag[-READ_LEN:])))

    # 4 host-derived pairs (removed by decontamination)
    for i in range(4):
        start = 1000 + i * 2000
        frag = host[start:start + INSERT]
        r1.append(fastq_record(f"host_pair_{i}/1", frag[:READ_LEN]))
        r2.append(fastq_record(f"host_pair_{i}/2", revcomp(frag[-READ_LEN:])))

    # 4 non-host singletons (survive)
    for i in range(4):
        s.append(fastq_record(f"nonhost_single_{i}", random_seq(rng, READ_LEN)))

    # 6 host-derived singletons (removed)
    for i in range(6):
        start = 20000 + i * 3000
        s.append(fastq_record(f"host_single_{i}", host[start:start + READ_LEN]))

    for name, records in [("R1", r1), ("R2", r2), ("S", s)]:
        out = OUT_DIR / f"singleton_test_{name}.fastq.gz"
        # mtime=0 keeps the gzip output byte-identical across regenerations
        with gzip.GzipFile(out, "wb", mtime=0) as handle:
            handle.write("".join(records).encode())
        print(f"wrote {out} ({len(records)} reads)")


if __name__ == "__main__":
    main()
