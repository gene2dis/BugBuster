#!/usr/bin/env python3
"""Split GTDB_TK_BATCH's pooled GTDB-Tk summaries back into per-sample tables.

GTDB_TK_BATCH classifies every sample's bins in one GTDB-Tk run, with each bin
copied as <sample>__<bin>.fa. This script writes <sample>_gtdbtk_<domain>.tsv
for every sample of the batch from a pooled gtdbtk.<domain>.summary.tsv, with
the sample prefix removed from user_genome.

Usage:
    gtdbtk_split.py --samples S1 S2 ... --domain bac120 --summary gtdbtk.bac120.summary.tsv
    gtdbtk_split.py --samples S1 S2 ... --domain ar53 [--summary <missing file>]

The sample id is the text before the first '__' (sample ids may not contain
'__': INPUT_CHECK rejects them). A bin name without '__', or whose prefix is not
a sample of the batch, is a naming bug and stops the script (exit 1) instead of
dropping the bin's row. Without a summary file (GTDB-Tk writes none for a
domain no genome was placed in) every sample gets a header-only table.
"""

import argparse
import os
import sys

# GTDB-Tk 2.8 summary columns, for the header-only tables of a domain with no
# summary file
SUMMARY_HEADER = "\t".join([
    "user_genome", "classification", "closest_genome_reference", "closest_genome_reference_radius",
    "closest_genome_taxonomy", "closest_genome_ani", "closest_genome_af", "closest_placement_reference",
    "closest_placement_radius", "closest_placement_taxonomy", "closest_placement_ani", "closest_placement_af",
    "pplacer_taxonomy", "classification_method", "note",
    "other_related_references(genome_id,species_name,radius,ANI,AF)", "msa_percent", "translation_table",
    "red_value", "warnings",
]) + "\n"


def split(summary, samples):
    """(header, {sample: [lines]}) of a pooled summary; exits 1 on a bin name
    whose sample prefix is missing or not in samples."""
    with open(summary) as fh:
        header = fh.readline()
        lines = [line for line in fh if line.strip()]
    rows = {sample: [] for sample in samples}
    for line in lines:
        name, sep, rest = line.partition("\t")
        sample, marker, bin_name = name.partition("__")
        if not marker or not bin_name:
            sys.exit(f"ERROR: {summary}: bin '{name}' has no '<sample>__' prefix; "
                     "GTDB_TK_BATCH names every pooled bin <sample>__<bin>")
        if sample not in rows:
            sys.exit(f"ERROR: {summary}: bin '{name}' has the prefix '{sample}', which is not a sample of "
                     f"this batch ({', '.join(samples)})")
        rows[sample].append(bin_name + sep + rest)
    return header, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", required=True, help="sample ids of the batch")
    parser.add_argument("--domain", required=True, choices=["bac120", "ar53"])
    parser.add_argument("--summary", help="pooled gtdbtk.<domain>.summary.tsv (missing or absent: header-only tables)")
    args = parser.parse_args()

    samples = list(dict.fromkeys(args.samples))
    if args.summary and os.path.exists(args.summary):
        header, rows = split(args.summary, samples)
    else:
        header, rows = SUMMARY_HEADER, {sample: [] for sample in samples}
    for sample in samples:
        with open(f"{sample}_gtdbtk_{args.domain}.tsv", "w") as out:
            out.write(header)
            out.writelines(rows[sample])
        print(f"{sample}_gtdbtk_{args.domain}.tsv: {len(rows[sample])} bins")


if __name__ == "__main__":
    main()
