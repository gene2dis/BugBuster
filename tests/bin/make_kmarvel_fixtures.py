#!/usr/bin/env python3
"""Generate tests/data/kmarvel/db_synthetic/ and the K-MARVEL test reads
(tests/modules/kmarvel.nf.test; roadmap Phase 25, RV25; review R8-06).

  db_synthetic/CARD_protein_all.fasta  three made-up 300-residue proteins in
                                       CARD's protein FASTA header layout
                                       (gb|<protein acc>|ARO:<id>|<name> [...])
  db_synthetic/metadata.tsv            their rows in CARD 4.0.1's metadata.tsv
                                       columns (the file K-MARVEL reads)
  db_synthetic/snps.txt                the snps.txt header only (no variant
                                       model)
  kmarvel_sim_R{1,2}.fastq.gz          200 error-free 150 bp pairs from a
                                       back-translated copy of SynR-1 with
                                       300 bp random flanks (~40x), plus 200
                                       random pairs

No DB_VERSION and no BLAST index: the module indexes a copy of the FASTA in
the task and records its md5 (the --custom_kmarvel_db path). The zero-ARG
case uses tests/data/reads/test_sample1.

Nothing here is CARD data (CARD's terms exclude commercial use; roadmap
decision RV8, synthetic CARD fixtures): names, accessions and sequences are
made up, only the layouts follow CARD 4.0.1. Byte-stable: fixed seed, gzip
mtime 0. Python standard library only; no network.
"""

import gzip
import random
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "data" / "kmarvel"
DB = OUT / "db_synthetic"
AA = "ACDEFGHIKLMNPQRSTVWY"
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
META_HEADER = ("ARO Accession\tCVTERM ID\tModel Sequence ID\tModel ID\tModel Name\tARO Name\t"
               "Protein Accession\tDNA Accession\tAMR Gene Family\tDrug Class\tResistance Mechanism\t"
               "CARD Short Name\n")
SNPS_HEADER = "Accession\tName\tModel Type\tParameter Type\tMutations\tCARD Short Name\tsource\tcitation\n"
COMP = str.maketrans("ACGT", "TGCA")


def codon_table():
    table, i = {}, 0
    for a in "TCAG":
        for b in "TCAG":
            for c in "TCAG":
                table.setdefault(CODE[i], []).append(a + b + c)
                i += 1
    return table


def revcomp(s):
    return s.translate(COMP)[::-1]


def write_gz(path, text):
    with open(path, "wb") as fh:
        with gzip.GzipFile(fileobj=fh, mode="wb", mtime=0) as gz:
            gz.write(text.encode())


def main():
    rng = random.Random(25)
    codons = codon_table()
    prots = ["M" + "".join(rng.choice(AA) for _ in range(299)) for _ in range(3)]
    DB.mkdir(parents=True, exist_ok=True)

    with open(DB / "CARD_protein_all.fasta", "w") as f:
        for n, p in enumerate(prots, 1):
            f.write(f">gb|SYN{n:05d}.1|ARO:900000{n}|SynR-{n} [synthetic organism]\n")
            f.writelines(p[k:k + 60] + "\n" for k in range(0, len(p), 60))
    with open(DB / "metadata.tsv", "w") as f:
        f.write(META_HEADER)
        for n in range(1, 4):
            f.write(f"ARO:900000{n}\t9900{n}\t990{n}\t99{n}\tsynthetic resistance protein SynR-{n}\t"
                    f"synthetic resistance protein SynR-{n}\tSYN{n:05d}.1\tSYNN{n:05d}.1\t"
                    f"synthetic resistance family {n}\tsynthetic antibiotic\tantibiotic inactivation\tSynR-{n}\n")
    (DB / "snps.txt").write_text(SNPS_HEADER)

    gene = "".join(rng.choice(codons[a]) for a in prots[0]) + "TAA"
    flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    genome = flank(300) + gene + flank(300)
    r1, r2 = [], []
    rec = lambda name, s: f"@{name}\n{s}\n+\n{'I' * len(s)}\n"
    for k in range(200):
        st = rng.randrange(0, len(genome) - 300)
        frag = genome[st:st + 300]
        if rng.random() < 0.5:
            frag = revcomp(frag)
        r1.append(rec(f"arg{k}/1", frag[:150]))
        r2.append(rec(f"arg{k}/2", revcomp(frag)[:150]))
    for k in range(200):
        frag = flank(300)
        r1.append(rec(f"rnd{k}/1", frag[:150]))
        r2.append(rec(f"rnd{k}/2", revcomp(frag)[:150]))
    write_gz(OUT / "kmarvel_sim_R1.fastq.gz", "".join(r1))
    write_gz(OUT / "kmarvel_sim_R2.fastq.gz", "".join(r2))


if __name__ == "__main__":
    main()
