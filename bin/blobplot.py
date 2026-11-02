#!/usr/bin/env python3
"""
Phylum blob plots: one figure per assembly (<assembly>_Phylum_blob_plot.pdf/.svg/.png + .csv), from the BlobTools
tables (*Blob_table*) in the working directory.

Left: each contig at its GC content and coverage (log scale), its area scaled by its length, coloured by phylum
(the 8 phyla with the most contig length, the rest as 'other'; 'no-hit' is grey and takes no colour). Right: each
phylum's share of the assembly's contig length. Coverage is BlobTools' read coverage (bam0; under co-assembly the
pooled reads of every sample). Drawn in the shared style of figure_style.py.
"""

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style as fs  # noqa: E402

fs.setup()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

TOP_PHYLA = 8


def read_blob_table(path):
    df = pd.read_csv(path, sep="\t", skiprows=10)
    return df[["# name", "length", "GC", "bam0", "phylum.t.9%s"]].rename(columns={
        "# name": "contig", "length": "length_bp", "GC": "gc", "bam0": "coverage", "phylum.t.9%s": "phylum"})


def label_phyla(data, top=TOP_PHYLA):
    """The top phyla by contig length keep their name; the rest are 'other'; 'no-hit' stays."""
    span = data[data["phylum"] != "no-hit"].groupby("phylum")["length_bp"].sum().sort_values(ascending=False,
                                                                                            kind="stable")
    keep = list(span.index[:top])
    data["group"] = np.where(data["phylum"].isin(keep) | (data["phylum"] == "no-hit"), data["phylum"], "other")
    order = keep + [g for g in ("other", "no-hit") if g in set(data["group"])]
    return data, order


def blob_figure(data, order, group_title, share_label, share):
    """Scatter of GC against coverage (area ~ length) and a bar of each group's share. share: {group: %}."""
    pal = fs.category_colors(order)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=fs.figsize("double", 80), layout="constrained",
                                 gridspec_kw={"width_ratios": [1.7, 1]})
    cov = data["coverage"].where(data["coverage"] > 0)
    floor = cov.min() / 2 if cov.notna().any() else 0.1
    y = cov.fillna(floor)
    size = 1 + 40 * data["length_bp"] / max(data["length_bp"].max(), 1)
    for g in reversed(order):  # the largest groups on top
        m = data["group"] == g
        if m.any():
            ax.scatter(data.loc[m, "gc"] * 100, y[m], s=size[m], color=pal[g], alpha=0.75, lw=0.2,
                       edgecolors="#555555" if g == "no-hit" else "none", label=g, rasterized=len(data) > 5000)
    ax.set_yscale("log")
    ax.set_xlim(max(0, data["gc"].min() * 100 - 5) if len(data) else 0, min(100, data["gc"].max() * 100 + 5)
                if len(data) else 100)
    ax.set_xlabel("GC content (%)")
    ax.set_ylabel("coverage (log scale)")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    yy = np.arange(len(order))
    bx.barh(yy, [share.get(g, 0) for g in order], color=[pal[g] for g in order], height=0.7, lw=0.3,
            edgecolor=["#555555" if g == "no-hit" else "none" for g in order])
    for i, g in enumerate(order):
        bx.text(share.get(g, 0) + 1, i, f"{share.get(g, 0):.1f}", va="center", fontsize=fs.SMALL_PT)
    bx.set_yticks(yy, order)
    bx.set_ylim(len(order) - 0.5, -0.5)
    bx.set_xlim(0, max(share.values(), default=1) * 1.18)
    bx.tick_params(axis="y", length=0)
    bx.set_xlabel(share_label)
    handles = [Line2D([], [], ls="", marker="o", markersize=4, color=pal[g], markeredgewidth=0.3,
                      markeredgecolor="#555555" if g == "no-hit" else pal[g], label=g) for g in order]
    fig.legend(handles=handles, title=group_title, loc="outside lower center", ncol=min(5, len(order)),
               handletextpad=0.3, columnspacing=1.0)
    return fig


def main():
    files = sorted(glob.glob("*Blob_table*"))
    if not files:
        print("ERROR: No blob table files found (*Blob_table*)", file=sys.stderr)
        sys.exit(1)
    for path in files:
        assembly = Path(path).name.split("_Blob_table")[0]
        data, order = label_phyla(read_blob_table(path))
        if data.empty:
            fs.placeholder(".", f"{assembly}_Phylum_blob_plot", f"{assembly}: no contigs in the BlobTools table")
            print(f"{assembly}_Phylum_blob_plot: no contigs")
            continue
        span = data.groupby("group")["length_bp"].sum()
        share = (span / span.sum() * 100).to_dict()
        fig = blob_figure(data, order, "phylum", "% of contig length", share)
        src = data[["contig", "length_bp", "gc", "coverage", "phylum", "group"]].assign(assembly=assembly)
        fs.save(fig, ".", f"{assembly}_Phylum_blob_plot", src)
        print(f"{assembly}_Phylum_blob_plot: {len(data)} contigs, phyla shown: {', '.join(order)}")


if __name__ == "__main__":
    main()
