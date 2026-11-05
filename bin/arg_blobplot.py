#!/usr/bin/env python3
"""
ARG blob plots: one figure per assembly (<assembly>_ARG_blob_plot.pdf/.svg/.png + .csv) of the ARG-carrying
contigs in Contig_tax_and_arg_prediction.tsv.

One point per contig (the table has one row per ARG hit): a contig whose hits fall in several ARG classes is
'multiple classes'. Colour: the 8 ARG classes with the most contigs, the rest as 'other'. Left: GC content against
coverage (log scale), area scaled by contig length. Right: each class's share of the assembly's ARG-carrying
contigs. An assembly without ARG hits gets a figure saying so. Drawn with blobplot.py's layout and the shared
style of figure_style.py.

Python port of ARG_blob_plot.R.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style as fs  # noqa: E402
from blobplot import blob_figure  # noqa: E402

TOP_CLASSES = 8
MULTI = "multiple classes"


def contig_classes(rows):
    """One row per ARG-carrying contig: its ARG class, or 'multiple classes'."""
    hits = rows.dropna(subset=["predicted_ARG_class"])
    hits = hits[hits["predicted_ARG_class"].astype(str).str.strip() != ""]
    if hits.empty:
        return pd.DataFrame(columns=["contig", "length_bp", "gc", "coverage", "arg_class", "arg_genes"])
    g = hits.groupby("contig_name", sort=False)
    out = g.agg(length_bp=("contig_length", "first"), gc=("GC", "first"), coverage=("coverage", "first"),
                arg_genes=("ARG_gene", lambda x: ";".join(sorted({str(v) for v in x if pd.notna(v)}))))
    classes = g["predicted_ARG_class"].agg(lambda x: sorted({str(v) for v in x}))
    out["arg_class"] = [c[0] if len(c) == 1 else MULTI for c in classes]
    out["arg_classes"] = [";".join(c) for c in classes]
    return out.rename_axis("contig").reset_index()


def main():
    try:
        table = pd.read_csv("Contig_tax_and_arg_prediction.tsv", sep="\t", dtype={"Sample": str}, low_memory=False)
    except FileNotFoundError:
        print("ERROR: Contig_tax_and_arg_prediction.tsv not found", file=sys.stderr)
        sys.exit(1)
    for assembly, rows in table.groupby("Sample", sort=True):
        key = f"{assembly}_ARG_blob_plot"
        data = contig_classes(rows)
        if data.empty:
            fs.placeholder(".", key, f"{assembly}: no ARG-carrying contigs")
            print(f"{key}: no ARG-carrying contigs")
            continue
        counts = data[data["arg_class"] != MULTI]["arg_class"].value_counts()
        keep = list(counts.sort_index(kind="stable").sort_values(ascending=False, kind="stable").index[:TOP_CLASSES])
        data["group"] = np.where(data["arg_class"].isin(keep) | (data["arg_class"] == MULTI), data["arg_class"],
                                 "other")
        order = keep + [g for g in (MULTI, "other") if g in set(data["group"])]
        share = (data["group"].value_counts() / len(data) * 100).to_dict()
        fig = blob_figure(data, order, "ARG class (DeepARG)", "% of ARG-carrying contigs", share)
        fs.save(fig, ".", key, data[["contig", "length_bp", "gc", "coverage", "arg_genes", "arg_classes", "group"]]
                .assign(assembly=assembly))
        print(f"{key}: {len(data)} ARG-carrying contigs, classes shown: {', '.join(order)}")


if __name__ == "__main__":
    main()
