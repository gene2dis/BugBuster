#!/usr/bin/env python3
"""Generate tests/data/arg_clustering/ (roadmap Phase 25, RV25), the inputs of
the bin-level ARG chain BIN_ARG_SPLIT -> ARG_FASTA_FORMATTER -> CLUSTERING in
tests/modules/arg_fasta_formatter.nf.test and clustering.nf.test.

Two samples, each in the layout DEEPARG_CONTIGS and BIN_CONTIG_LABELS give
BIN_ARG_SPLIT:
  <s>_contig_bins.tsv           contig -> bin
  <s>_proteins.faa              Pyrodigal-style proteins of the assembly
                                (<contig>_<n> ids, '#'-separated header)
  <s>_deep_arg.out.mapping.ARG  DeepARG's header and its hits

Proteins (made up, seeded random 300-mers):
  P1   in s1 bin1 (k141_1_1) and, identical, in s2 bin1 (k141_7_1): one
       cluster at every tier, 100 % included
  P1v  P1 with 21 of its 300 residues changed (93 % identity), s2 bin2
       (k141_9_1): with P1 at 90 %, apart at 95 % and above
  P2   unrelated, s1 bin1 (k141_1_2): a cluster of its own at every tier
  P3   s1 bin1 (k141_2_1), no DeepARG hit: not written by the formatter
s1 bin2 has a protein but no hit, so the formatter writes no file for it.

Nothing here is DeepARG or CARD data: names, classes and sequences are made
up. Deterministic: re-running reproduces byte-identical files.
"""

import random
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "arg_clustering"
AA = "ACDEFGHIKLMNPQRSTVWY"
HEADER = ("#ARG\tquery-start\tquery-end\tread_id\tpredicted_ARG-class\tbest-hit\tprobability"
          "\tidentity\talignment-length\talignment-bitscore\talignment-evalue\tcounts\n")


def protein(rng, n=300):
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def mutate(rng, seq, n):
    seq = list(seq)
    for pos in rng.sample(range(1, len(seq)), n):
        seq[pos] = rng.choice([a for a in AA if a != seq[pos]])
    return "".join(seq)


def faa(records):
    out = []
    for i, (gene, seq) in enumerate(records, 1):
        out.append(f">{gene} # 1 # {3 * len(seq) + 3} # 1 # ID={i}_1;partial=00;start_type=ATG\n")
        out.extend(seq[k:k + 60] + "\n" for k in range(0, len(seq), 60))
    return "".join(out)


def hit(gene, arg, arg_class, seq):
    n = len(seq)
    return (f"{arg}\t1\t{n}\t{gene}\t{arg_class}\t{arg}|FEATURES|{arg.lower()}|{arg_class}|{arg.lower()}"
            f"\t0.99\t100.0\t{n}\t600.0\t1e-150\t1\n")


def main():
    rng = random.Random(25)
    p1 = protein(rng)
    p2 = protein(rng)
    p3 = protein(rng)
    p4 = protein(rng)
    p1v = mutate(rng, p1, 21)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    samples = {
        "s1": {
            "bins": [("k141_1", "bin1"), ("k141_2", "bin1"), ("k141_3", "bin2")],
            "proteins": [("k141_1_1", p1), ("k141_1_2", p2), ("k141_2_1", p3), ("k141_3_1", p4)],
            "hits": [("k141_1_1", "SYNA", "beta-lactam", p1), ("k141_1_2", "SYNB", "tetracycline", p2)],
        },
        "s2": {
            "bins": [("k141_7", "bin1"), ("k141_9", "bin2")],
            "proteins": [("k141_7_1", p1), ("k141_9_1", p1v)],
            "hits": [("k141_7_1", "SYNA", "beta-lactam", p1), ("k141_9_1", "SYNA", "beta-lactam", p1v)],
        },
    }
    for s, d in samples.items():
        (OUT_DIR / f"{s}_contig_bins.tsv").write_text(
            "contig\tbin\n" + "".join(f"{c}\t{b}\n" for c, b in d["bins"]))
        (OUT_DIR / f"{s}_proteins.faa").write_text(faa(d["proteins"]))
        (OUT_DIR / f"{s}_deep_arg.out.mapping.ARG").write_text(
            HEADER + "".join(hit(*h) for h in d["hits"]))


if __name__ == "__main__":
    main()
