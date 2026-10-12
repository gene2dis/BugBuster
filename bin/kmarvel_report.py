#!/usr/bin/env python3
"""
K-MARVEL summary report.

Usage: kmarvel_report.py <sample>_kmarvel.tsv [...]

Merges the per-sample K-MARVEL tables (one row per assembled ARG, written by
the KMARVEL module) into:
  KMARVEL_summary_long.tsv             one row per sample and ARG, key columns
  KMARVEL_supporting_reads_matrix.tsv  ARG (CARD ARO) x sample, supporting read
                                       pairs (0 = not found)
  KMARVEL_gene_family_matrix.tsv       AMR gene family x sample, supporting read
                                       pairs summed over the family's ARGs
  KMARVEL_drug_class_matrix.tsv        drug class x sample, number of ARGs found
                                       (an ARG counts once for each ';'-separated
                                       class it confers resistance to)
K-MARVEL's supporting_reads counts read pairs: both mates of a fragment count
once (checked on simulated pairs with /1 /2, Illumina "1:N:0" and identical
mate headers). The allele it reports among close variants of one gene (TEM-1
reads called TEM-117, TEM-192 or TEM-246 depending on the read subset) is not
stable, so compare samples at gene-family level. A gene family is kept exactly
as CARD lists it: a combined family such as "AAC(6');AAC(6')-Ib-cr" is one row,
so no read is counted twice.
A sample with no ARGs (header-only table) is kept as an all-zero column.
Counts are not normalised: K-MARVEL can run without ARGs-OAP, which provides
the cell counts the KARGA/KARGVA report uses.
"""

import sys
from pathlib import Path

import pandas as pd

SUFFIX = "_kmarvel.tsv"
LONG_COLUMNS = [
    "aro_id", "aro_full_name", "gene_family", "drug_class",
    "resistance_mechanism", "model", "identity", "assembly_coverage",
    "supporting_reads", "median_kmer_depth", "blast_predicted_aro_name",
    "blast_e_value", "blast_bit_score", "mutation_data",
]


def sample_name(path):
    name = Path(path).name
    if not name.endswith(SUFFIX):
        sys.exit(f"ERROR: {name} does not end with {SUFFIX}")
    return name[: -len(SUFFIX)]


def main(paths):
    if not paths:
        sys.exit("ERROR: no K-MARVEL tables given")

    samples = []
    tables = []
    for path in paths:
        sample = sample_name(path)
        samples.append(sample)
        # Check the header first, then read only the columns used: the
        # read_data, sequence and CIGAR columns grow with read depth
        header = pd.read_csv(path, sep="\t", nrows=0).columns
        missing = [c for c in LONG_COLUMNS if c not in header]
        if missing:
            sys.exit(f"ERROR: {path} lacks columns: {', '.join(missing)}")
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                         usecols=LONG_COLUMNS)[LONG_COLUMNS].copy()
        df.insert(0, "sample", sample)
        tables.append(df)

    samples = sorted(set(samples))
    long = pd.concat(tables, ignore_index=True)
    long["supporting_reads"] = pd.to_numeric(long["supporting_reads"])
    long = long.sort_values(["sample", "aro_id"]).reset_index(drop=True)
    long.to_csv("KMARVEL_summary_long.tsv", sep="\t", index=False)

    reads = (
        long.pivot_table(index=["aro_id", "aro_full_name"], columns="sample",
                         values="supporting_reads", aggfunc="sum", fill_value=0)
        .reindex(columns=samples, fill_value=0)
    )
    reads.columns.name = None
    reads.to_csv("KMARVEL_supporting_reads_matrix.tsv", sep="\t")

    family = (
        long.pivot_table(index="gene_family", columns="sample",
                         values="supporting_reads", aggfunc="sum", fill_value=0)
        .reindex(columns=samples, fill_value=0)
    )
    family.columns.name = None
    family.to_csv("KMARVEL_gene_family_matrix.tsv", sep="\t")

    classes = long[["sample", "aro_id", "drug_class"]].copy()
    classes["drug_class"] = classes["drug_class"].str.split(";")
    classes = classes.explode("drug_class")
    classes["drug_class"] = classes["drug_class"].str.strip()
    classes = classes[classes["drug_class"] != ""]
    drug = (
        classes.groupby(["drug_class", "sample"])["aro_id"].nunique()
        .unstack(fill_value=0)
        .reindex(columns=samples, fill_value=0)
    )
    drug.index.name = "drug_class"
    drug.columns.name = None
    drug.to_csv("KMARVEL_drug_class_matrix.tsv", sep="\t")

    print(f"{len(samples)} samples, {len(long)} ARG calls, "
          f"{len(reads)} distinct ARGs, {len(family)} gene families, "
          f"{len(drug)} drug classes")


if __name__ == "__main__":
    main(sys.argv[1:])
