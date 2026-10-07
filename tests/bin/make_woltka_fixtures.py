#!/usr/bin/env python3
"""Generate the tiny deterministic WoLr2-layout fixture used by the Woltka
read-branch tests (design doc T8a, Section 4.6.2).

Writes tests/data/woltka/ (sequences, coordinates, maps, reads). The Bowtie2
index is built afterwards by tests/bin/make_woltka_fixtures.sh, which calls
this script first.

Database layout mirrors https://ftp.microbio.me/pub/wol2/ verbatim (only the
files the backend reads), so the fixture doubles as a --custom_woltka_db:

  db/databases/bowtie2/WoLr2.*.bt2l     built by the .sh wrapper
  db/proteins/coords.txt.xz (+ .md5 of the UNCOMPRESSED content, as on the FTP)
  db/proteins/length.map.xz (+ .md5)
  db/function/kegg/{orf-to-ko.map.xz(+.md5),ko-to-ec.map,ko-to-cog.map,ko_name.txt}
  db/function/metacyc/{orf-to-protein.map.xz(+.md5),protein-to-enzrxn.map,
      enzrxn-to-reaction.map,reaction-to-pathway.map,pathway_name.txt}
  db/function/pfam/{orf-to-pfam.map.xz(+.md5),pfam_name.txt}
  db/DB_VERSION
  db/genomes.fna                         index-build input (not part of WoLr2)

Two 4 kb genomes. G000000002's first 1500 bp is an exact copy of
G000000001's, so reads drawn there align equally well to both genomes:
  G000000001_1  101-1400 (+)  shared region   KOs K00001 K00002, Pfam dup row
  G000000001_2 1601-2800 (-)  unique          K00003
  G000000001_3 3001-3900 (+)  unique          unannotated in every map
  G000000002_1  101-1400 (+)  shared region   K00001
  G000000002_2 2001-3500 (+)  unique          no KO; MetaCyc + Pfam only

Map rows deliberately exercise the term-set composition (woltka collapse
inflates both, verified 2026-10-07 on 0.1.7):
  - Pfam: G000000001_1 lists PF00002.1 three times (WoLr2 repeats domains)
  - KEGG: G000000001_1 -> K00001 + K00002, both -> EC 1.1.1.1 (EC once/gene)
  - MetaCyc: G000000001_1 -> P1 + P2; P1 -> E1 + E2 -> both R1; P2 -> E3 -> R2;
    R1 -> PWY-1; R2 -> PWY-1 + PWY-2 (PWY-1 must count once per gene)

Reads (exact 100-mers, Q40, ~300 bp inserts, proper FR pairs):
  sample  'wtest': 4 pairs inside G000000001_2 (unique), 3 pairs inside the
          shared G*_1 region (ambiguous between the genomes), 2 pairs inside
          G000000002_2, 1 pair in the G000000001 2801-3000 intergenic gap,
          plus 2 singletons inside G000000001_3.

Also writes tests/data/woltka/wtest.ags.tsv (GE 2.0) for the aggregation
tests; the .sh wrapper adds the alignment and the WOLTKA_CLASSIFY outputs.
"""

import gzip
import hashlib
import lzma
import random
from pathlib import Path

READ_LEN = 100
INSERT = 300
QUAL = "I" * READ_LEN

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "data" / "woltka"
DB = OUT / "db"

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def random_seq(rng, length):
    while True:
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        if "G" * 10 not in seq and "C" * 10 not in seq:
            return seq


def write_xz_with_md5(path, text):
    """Write <path>.xz plus <path>.md5 holding the md5 of the uncompressed
    text, which is how the WoLr2 FTP checksums its .xz files."""
    path.parent.mkdir(parents=True, exist_ok=True)
    data = text.encode()
    with lzma.open(f"{path}.xz", "wb", preset=6) as fh:
        fh.write(data)
    Path(f"{path}.md5").write_text(hashlib.md5(data).hexdigest() + "\n")


def write_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def pair_from(rng, genome, start, end, name):
    """A proper FR pair whose 300 bp fragment lies inside [start, end] (1-based,
    inclusive), so both mates fall >= 80 % inside the region."""
    frag_start = rng.randint(start - 1, end - INSERT)  # 0-based
    frag = genome[frag_start:frag_start + INSERT]
    r1 = frag[:READ_LEN]
    r2 = revcomp(frag[-READ_LEN:])
    return (f"@{name}/1\n{r1}\n+\n{QUAL}\n", f"@{name}/2\n{r2}\n+\n{QUAL}\n")


def single_from(rng, genome, start, end, name):
    s = rng.randint(start - 1, end - READ_LEN)
    return f"@{name}\n{genome[s:s + READ_LEN]}\n+\n{QUAL}\n"


def main():
    rng = random.Random(8)
    g1 = random_seq(rng, 4000)
    g2 = g1[:1500] + random_seq(rng, 2500)

    write_text(DB / "genomes.fna", f">G000000001\n{g1}\n>G000000002\n{g2}\n")

    # WoLr2 coords.txt: '>genome' then 'idx<TAB>start<TAB>end' (start > end
    # marks the minus strand), lengths in length.map are end - start + 1
    orfs = {
        "G000000001": [(1, 101, 1400), (2, 2800, 1601), (3, 3001, 3900)],
        "G000000002": [(1, 101, 1400), (2, 2001, 3500)],
    }
    coords, lengths = [], []
    for genome, rows in orfs.items():
        coords.append(f">{genome}")
        for idx, a, b in rows:
            coords.append(f"{idx}\t{a}\t{b}")
            lengths.append(f"{genome}_{idx}\t{abs(b - a) + 1}")
    write_xz_with_md5(DB / "proteins" / "coords.txt", "\n".join(coords) + "\n")
    write_xz_with_md5(DB / "proteins" / "length.map", "\n".join(lengths) + "\n")

    kegg = DB / "function" / "kegg"
    write_xz_with_md5(kegg / "orf-to-ko.map",
                      "G000000001_1\tK00001\tK00002\n"
                      "G000000001_2\tK00003\n"
                      "G000000002_1\tK00001\n")
    write_text(kegg / "ko-to-ec.map",
               "K00001\t1.1.1.1\n"
               "K00002\t1.1.1.1\t2.7.7.7\n"
               "K00003\t3.1.1.1\n")
    write_text(kegg / "ko-to-cog.map",
               "K00001\tCOG0001\n"
               "K00002\tCOG0002\n"
               "K00003\tCOG0001\tCOG0003\n")
    write_text(kegg / "ko_name.txt",
               "K00001\tE1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]\n"
               "K00002\tfixture enzyme two\n"
               "K00003\tfixture enzyme three\n")

    metacyc = DB / "function" / "metacyc"
    write_xz_with_md5(metacyc / "orf-to-protein.map",
                      "G000000001_1\tP1\tP2\n"
                      "G000000002_2\tP3\n")
    write_text(metacyc / "protein-to-enzrxn.map",
               "P1\tE1\tE2\n"
               "P2\tE3\n"
               "P3\tE3\n")
    write_text(metacyc / "enzrxn-to-reaction.map",
               "E1\tR1\n"
               "E2\tR1\n"
               "E3\tR2\n")
    write_text(metacyc / "reaction-to-pathway.map",
               "R1\tPWY-1\n"
               "R2\tPWY-1\tPWY-2\n")
    write_text(metacyc / "pathway_name.txt",
               "PWY-1\tfixture pathway one\n"
               "PWY-2\tfixture pathway two (<i>E. coli</i>)\n")

    pfam = DB / "function" / "pfam"
    write_xz_with_md5(pfam / "orf-to-pfam.map",
                      "G000000001_1\tPF00001.1\tPF00002.1\tPF00002.1\tPF00002.1\n"
                      "G000000001_2\tPF00001.1\n"
                      "G000000002_2\tPF00003.1\n")
    write_text(pfam / "pfam_name.txt",
               "PF00001.1\t7tm_1\n"
               "PF00002.1\t7tm_2\n"
               "PF00003.1\t7tm_3\n")

    write_text(DB / "DB_VERSION", "WoLr2 test fixture (tests/bin/make_woltka_fixtures.py)\n")

    r1, r2, s = [], [], []
    n = 0
    for genome, start, end, count in [
        (g1, 1601, 2800, 4),   # G000000001_2, unique
        (g1, 101, 1400, 3),    # shared region: G000000001_1 / G000000002_1
        (g2, 2001, 3500, 2),   # G000000002_2, unique
    ]:
        for _ in range(count):
            n += 1
            a, b = pair_from(rng, genome, start, end, f"read{n}")
            r1.append(a)
            r2.append(b)
    # one pair spanning the G000000001 intergenic gap 2801-3000: each mate
    # overlaps no ORF by >= 80 % -> aligned but unassigned
    n += 1
    frag = g1[2750:2750 + INSERT]
    r1.append(f"@read{n}/1\n{frag[:READ_LEN]}\n+\n{QUAL}\n")
    r2.append(f"@read{n}/2\n{revcomp(frag[-READ_LEN:])}\n+\n{QUAL}\n")
    for i in range(2):
        s.append(single_from(rng, g1, 3001, 3900, f"single{i + 1}"))

    # MICROBECENSUS-shaped AGS table for the aggregation tests: GE = 2.0
    # makes CPGE = RPK / 2 hand-checkable
    write_text(OUT / "wtest.ags.tsv",
               "sample_id\taverage_genome_size_bp\tgenome_equivalents\ttotal_bases\n"
               "wtest\t3000000\t2.0\t6000000\n")

    reads = OUT / "reads"
    reads.mkdir(parents=True, exist_ok=True)
    for name, records in [("wtest_R1", r1), ("wtest_R2", r2), ("wtest_S", s)]:
        with gzip.GzipFile(reads / f"{name}.fastq.gz", "wb", mtime=0) as fh:
            fh.write("".join(records).encode())


if __name__ == "__main__":
    main()
