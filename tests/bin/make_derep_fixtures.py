#!/usr/bin/env python3
"""Generate tests/data/derep/, the fixtures of tests/modules/galah_cluster.nf.test and
tests/modules/coverm_genome.nf.test (--dereplicate_mags, roadmap Phase 18b, B5).

Two assemblies, 'A' and 'site_1' (an id with '_', which the cluster table must not split on):
- A_binette_bins/binette_bin1.fa: genome X (60 kb, two contigs); binette_bin2.fa: genome Y (60 kb);
- site_1_binette_bins/binette_bin1.fa: X with 1 % of its bases changed (the same species at 95 % ANI).
Binette-style CheckM2 reports: A bin1 90 % / 1 % (score 85 by completeness - 5 x contamination),
site_1 bin1 95 % / 3 % (score 80, but the more complete one): galah's completeness-5contamination formula
picks A's bin as X's representative. reads_R{1,2}.fastq.gz: 2,000 read pairs (2 x 100 bp, 300 bp inserts,
no errors) from Y only. Seeded: the same files on every run.
"""
import gzip
import random
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "data" / "derep"
rng = random.Random(18)
BASES = "ACGT"


def genome(n):
    return "".join(rng.choice(BASES) for _ in range(n))


def mutate(seq, rate):
    s = list(seq)
    for i in rng.sample(range(len(s)), int(len(s) * rate)):
        s[i] = rng.choice([b for b in BASES if b != s[i]])
    return "".join(s)


def write_fasta(path, contigs):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in contigs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


HEADER = ("Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\tCoding_Density\t"
          "Contig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tTotal_Contigs\t"
          "Max_Contig_Length\tAdditional_Notes\n")


def report_row(name, compl, cont, size):
    return f"{name}\t{compl}\t{cont}\tNeural Network (Specific Model)\t11\t0.87\t30000\t300.0\t{size}\t0.5\t60\t2\t30000\t\n"


x = genome(60000)
y = genome(60000)
xm = mutate(x, 0.01)
write_fasta(OUT / "A_binette_bins" / "binette_bin1.fa", [("A_k141_1", x[:30000]), ("A_k141_2", x[30000:])])
write_fasta(OUT / "A_binette_bins" / "binette_bin2.fa", [("A_k141_3", y[:30000]), ("A_k141_4", y[30000:])])
write_fasta(OUT / "site_1_binette_bins" / "binette_bin1.fa", [("S_k141_1", xm[:30000]), ("S_k141_2", xm[30000:])])
(OUT / "A_binette_quality_report.tsv").write_text(
    HEADER + report_row("binette_bin1", 90.0, 1.0, 60000) + report_row("binette_bin2", 99.0, 0.5, 60000))
(OUT / "site_1_binette_quality_report.tsv").write_text(HEADER + report_row("binette_bin1", 95.0, 3.0, 60000))

with gzip.GzipFile(OUT / "reads_R1.fastq.gz", "wb", mtime=0) as r1, \
        gzip.GzipFile(OUT / "reads_R2.fastq.gz", "wb", mtime=0) as r2:
    for i in range(2000):
        start = rng.randrange(0, len(y) - 300)
        frag = y[start:start + 300]
        q = "I" * 100
        r1.write(f"@r{i}/1\n{frag[:100]}\n+\n{q}\n".encode())
        r2.write(f"@r{i}/2\n{revcomp(frag[-100:])}\n+\n{q}\n".encode())
