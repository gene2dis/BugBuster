#!/usr/bin/env python3
"""Deterministic reads for the SUPER-FOCUS read-branch fixture (design doc T8b).

Source data (committed, not generated):
  tests/data/superfocus/source/proteins.faa   five real SUPER-FOCUS DB_90
      cluster proteins (diamond getseq of the figshare CC0 database, one per
      SEED subsystem PK 150/151/615/678/134) plus a synthetic copy of the
      PK 150 fumarate hydratase sequence under PK 151, so its reads tie
      between two subsystems (the 1/k split)
  tests/data/superfocus/sf_db/db/database_PKs.txt   the matching rows of the v1.8
      database_PKs.txt (two level-2 '-' subsystems under different level-1
      categories: the path-qualified accession case)

Writes tests/data/superfocus/reads/sftest_{R1,R2,S}.fastq.gz: each read is
a 150 nt window of a protein's reverse translation (R1 forward, its mate R2
the reverse complement of the next window; mates share the read id, as in
real Illumina output), plus unrelated pseudo-random reads that hit nothing.

Expected SUPER-FOCUS counts (-n 1), all hand-checkable:
  FumHan  3 pairs  -> 6 reads, tie: 3 to TCA Cycle (150), 3 to TCA cycle in plants (151)
  SDH1    2 pairs  -> 4 reads to TCA cycle in plants (151)
  CMSA    1 pair + 1 singleton -> 3 reads to Creatine and Creatinine Degradation (615)
  acoR    2 pairs  -> 4 reads to Dehydrogenase complexes (678)
  TM0313  1 pair   -> 2 reads to Sugar utilization in Thermotogales (134)
  random  1 pair + 1 singleton -> 3 reads, no hit
  => 22 input reads, 19 with an accepted hit.

Usage: python3 tests/bin/make_superfocus_fixtures.py   (stdlib only)
"""

import gzip
from pathlib import Path

FIX = Path(__file__).resolve().parents[2] / 'tests' / 'data' / 'superfocus'
READ_LEN = 150

# one fixed codon per amino acid (E. coli-preferred)
CODON = {
    'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC', 'Q': 'CAG', 'E': 'GAA',
    'G': 'GGC', 'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT',
    'P': 'CCG', 'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
}
COMPLEMENT = str.maketrans('ACGT', 'TGCA')

# (function tag in the header, number of read pairs, singleton windows)
PLAN = [('__150__TCA_Cycle__FumHan__', 3, 0),
        ('__151__TCA_cycle_in_plants__SDH1__', 2, 0),
        ('__615__Creatine_and_Creatinine_Degradation__CMSA__', 1, 1),
        ('__678__Dehydrogenase_complexes__acoR__', 2, 0),
        ('__134__Sugar_utilization_in_Thermotogales__TM0313__', 1, 0)]


def read_fasta(path):
    records, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:]
            records[name] = ''
        else:
            records[name] += line.strip()
    return records


def random_dna(n, state):
    """Linear congruential generator: reproducible on every platform."""
    out = []
    for _ in range(n):
        state = (state * 1103515245 + 12345) % 2**31
        out.append('ACGT'[(state >> 16) % 4])
    return ''.join(out), state


def fastq(name, seq):
    return f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n"


def write_gz(path, text):
    # mtime=0 keeps the file byte-stable across regenerations
    with open(path, 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0,
                                                 filename='') as handle:
        handle.write(text.encode())


def main():
    proteins = read_fasta(FIX / 'source' / 'proteins.faa')
    r1, r2, single = [], [], []
    fragment = 0
    for tag, pairs, singletons in PLAN:
        name = next(h for h in proteins if tag in h and not h.startswith('fixture_'))
        dna = ''.join(CODON[aa] for aa in proteins[name])
        for j in range(pairs):
            start = j * 2 * READ_LEN
            assert start + 2 * READ_LEN <= len(dna), (tag, j)
            fragment += 1
            r1.append(fastq(f"sf_frag{fragment} 1:N:0", dna[start:start + READ_LEN]))
            mate = dna[start + READ_LEN:start + 2 * READ_LEN].translate(COMPLEMENT)[::-1]
            r2.append(fastq(f"sf_frag{fragment} 2:N:0", mate))
        for k in range(singletons):
            start = len(dna) - READ_LEN - k * READ_LEN
            single.append(fastq(f"sf_single_{tag.split('__')[3]}{k + 1}", dna[start:start + READ_LEN]))

    state = 20261007
    fragment += 1
    seq, state = random_dna(READ_LEN, state)
    r1.append(fastq(f"sf_frag{fragment} 1:N:0", seq))
    seq, state = random_dna(READ_LEN, state)
    r2.append(fastq(f"sf_frag{fragment} 2:N:0", seq))
    seq, state = random_dna(READ_LEN, state)
    single.append(fastq("sf_single_random1", seq))

    reads = FIX / 'reads'
    reads.mkdir(parents=True, exist_ok=True)
    write_gz(reads / 'sftest_R1.fastq.gz', ''.join(r1))
    write_gz(reads / 'sftest_R2.fastq.gz', ''.join(r2))
    write_gz(reads / 'sftest_S.fastq.gz', ''.join(single))
    print(f"wrote {len(r1)} pairs + {len(single)} singletons to {reads}")


if __name__ == '__main__':
    main()
