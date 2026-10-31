"""Figures of bin/results_summary.py (PNG, PDF and SVG through FigureWriter)."""

import math

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from summary_common import NA_NOT_RUN, strip_gtdb  # noqa: E402
from summary_taxonomy import genus_rank, species_rank  # noqa: E402


def fontsize_for(n, lo=4, hi=9, scale=320):
    return max(lo, min(hi, scale / max(n, 1)))


def group_spans(cols, groups):
    spans, start = [], 0
    for i in range(1, len(cols) + 1):
        if i == len(cols) or groups[i] != groups[start]:
            spans.append((groups[start], start, i - 1))
            start = i
    return spans


class FigureWriter:
    def __init__(self, outdir):
        self.dir = outdir / "figures"
        self.dir.mkdir(parents=True, exist_ok=True)
        self.records = []

    def save(self, fig, key, title, caption):
        for ext in ("png", "pdf", "svg"):
            kw = {"dpi": 200} if ext == "png" else {}
            fig.savefig(self.dir / f"{key}.{ext}", bbox_inches="tight", facecolor="white", **kw)
        plt.close(fig)
        self.records.append({"key": key, "title": title, "caption": caption})


def label_levels(cols, groups, fontsize, inch_per_col):
    """Text row for each group name above the columns: the lowest row where it does not overlap a name already
    placed there (names wider than their group spill over the neighbours). Returns [(name, centre, row)]."""
    placed, out = [], []
    for g, start, end in group_spans(cols, groups):
        centre = (start + end) / 2
        half = len(str(g)) * fontsize * 0.62 / 72 / inch_per_col / 2 + 0.4
        lo, hi = centre - half, centre + half
        row = 0
        while any(r == row and not (hi < a or lo > b) for a, b, r in placed):
            row += 1
        placed.append((lo, hi, row))
        out.append((g, centre, row))
    return out


def place_group_labels(fig, ax, gax, cols, groups, fontsize):
    """Group names in gax (x shared with ax), each on the lowest text row where its measured extent does not
    overlap a name already there. Returns the number of rows used."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px_per_col = ax.get_window_extent(renderer).width / max(len(cols), 1)
    placed, texts = [], []
    for g, start, end in group_spans(cols, groups):
        centre = (start + end) / 2
        t = gax.text(centre, 0.05, str(g), ha="center", va="bottom", fontsize=fontsize, clip_on=False)
        # names must not move the layout: a name wider than its group would shrink the heatmap and void
        # the measured column width
        t.set_in_layout(False)
        half = t.get_window_extent(renderer).width / px_per_col / 2 + 0.5
        row = 0
        while any(r == row and not (centre + half < lo or centre - half > hi) for lo, hi, r in placed):
            row += 1
        placed.append((centre - half, centre + half, row))
        texts.append((t, row))
    n_levels = max(r for _, r in texts) + 1
    for t, row in texts:
        t.set_y(row + 0.05)
    gax.set_ylim(0, n_levels)
    return n_levels


def label_groups(ax, cols, groups, y0, step, fontsize, inch_per_col):
    """Group names above the columns (data coordinates: row k at y0 + k * step), staggered so none overlap."""
    for g, centre, row in label_levels(cols, groups, fontsize, inch_per_col):
        ax.text(centre, y0 + row * step, str(g), ha="center", va="bottom", fontsize=fontsize, clip_on=False)


def style_sample_ticks(labels_list, low, axis_labels):
    for lab in axis_labels:
        if lab.get_text() in low:
            lab.set_color("#c0392b")


def heatmap(fw, key, title, caption, mat, row_labels, cols, groups, low, hatch_cols, cbar_label,
            vmin=None, vmax=None, cmap="viridis", row_groups=None, row_group_title=""):
    """log10 heatmap, samples as columns. Group names sit in their own strip above the columns and the title
    is a figure title, so neither can overlap the other. row_groups (one label per row, e.g. phylum) adds a
    colour strip left of the rows with a legend below."""
    data = mat.reindex(columns=cols).astype(float).values
    with np.errstate(divide="ignore", invalid="ignore"):
        logd = np.where(data > 0, np.log10(data), np.nan)
    if np.all(np.isnan(logd)):
        return
    nrow, ncol = logd.shape
    row_h = 0.2 if nrow <= 60 else 0.15
    fs_x, fs_y = fontsize_for(ncol), max(4, min(9, row_h * 72 * 0.62))
    width = min(60, max(7, 0.2 * ncol + 4))
    heat_h = max(2.0, row_h * nrow)
    has_groups = any(groups)
    group_h = 0.5 if has_groups else 0.0
    legend_h = 0.0
    if row_groups is not None:
        legend_h = 0.3 + 0.22 * math.ceil(len(set(row_groups)) / 4)
    fig = plt.figure(figsize=(width, heat_h + group_h + legend_h + 2.6), layout="constrained")
    ncols_grid = 2 if row_groups is not None else 1
    gs = fig.add_gridspec(2 if has_groups else 1, ncols_grid,
                          height_ratios=[group_h, heat_h] if has_groups else [heat_h],
                          width_ratios=[0.18, width] if row_groups is not None else [1])
    ax = fig.add_subplot(gs[-1, -1])
    cm = plt.get_cmap(cmap)
    # drop the palest end so a low value never looks like the grey 'not detected' cell
    cm = matplotlib.colors.ListedColormap(cm(np.linspace(0.12, 1, 256)) if cmap.endswith("_r") else cm(np.linspace(0, 0.92, 256)))
    cm.set_bad("#f0f0f0")
    vmin = np.nanmin(logd) if vmin is None else vmin
    vmax = np.nanmax(logd) if vmax is None else vmax
    im = ax.imshow(np.ma.masked_invalid(logd), aspect="auto", cmap=cm, vmin=vmin, vmax=vmax,
                   interpolation="nearest")
    for j, s in enumerate(cols):
        if s in hatch_cols:
            ax.add_patch(Rectangle((j - 0.5, -0.5), 1, nrow, fill=False, hatch="///",
                                   edgecolor="#999999", linewidth=0))
    for _, start, end in group_spans(cols, groups)[:-1]:
        ax.axvline(end + 0.5, color="white", linewidth=2.5)
    ax.set_xticks(range(ncol))
    ax.set_xticklabels(cols, rotation=90, fontsize=fs_x)
    style_sample_ticks(cols, low, ax.get_xticklabels())
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    if row_groups is not None:
        sax = fig.add_subplot(gs[-1, 0], sharey=ax)
        names = list(dict.fromkeys(row_groups))
        palette = category_colors(names)
        sax.imshow(np.array([[names.index(g)] for g in row_groups]), aspect="auto", interpolation="nearest",
                   cmap=matplotlib.colors.ListedColormap([palette[n] for n in names]), vmin=-0.5,
                   vmax=len(names) - 0.5)
        sax.set_xticks([])
        sax.set_yticks(range(nrow))
        sax.set_yticklabels(row_labels, fontsize=fs_y)
        sax.tick_params(length=0)
        for sp in sax.spines.values():
            sp.set_visible(False)
        ax.tick_params(labelleft=False)
        handles = [Rectangle((0, 0), 1, 1, color=palette[n]) for n in names]
        fig.legend(handles, display_names(names), loc="outside lower center", ncol=4, fontsize=8,
                   frameon=False, title=row_group_title or None, title_fontsize=8)
    else:
        ax.set_yticks(range(nrow))
        ax.set_yticklabels(row_labels, fontsize=fs_y)
    if has_groups:
        gax = fig.add_subplot(gs[0, -1], sharex=ax)
        gax.axis("off")
        n_levels = place_group_labels(fig, ax, gax, cols, groups, fs_x + 1)
        gs.set_height_ratios([0.08 + 0.22 * n_levels, heat_h])
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, aspect=max(20, min(60, nrow / 2)))
    cb.set_label(cbar_label, fontsize=9)
    fig.suptitle(title, fontsize=12)
    fw.save(fig, key, title, caption)


def short_count(n):
    """10000000 -> '10M', 1500000 -> '1.5M', 250000 -> '250k'."""
    for div, suffix in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if n >= div:
            return f"{n / div:.3g}{suffix}"
    return f"{n:.0f}"


def category_colors(names):
    """A colour per category name, in the order given: tab10 for up to 10, else tab20 with its dark shades first
    (tab20 pairs each colour with a light shade of itself, so neighbours would look alike), then tab20b."""
    if len(names) <= 10:
        base = list(plt.get_cmap("tab10").colors)
    else:
        t20 = list(plt.get_cmap("tab20").colors)
        base = t20[0::2] + t20[1::2] + list(plt.get_cmap("tab20b").colors)
    return {n: base[i % len(base)] for i, n in enumerate(names)}


def fig_qc(fw, qc, cols, groups, low, min_pairs, depth_marks=()):
    if "raw_pairs" not in qc or qc["raw_pairs"].isna().all():
        return
    q = qc.reindex(cols)
    n = len(cols)
    fs = fontsize_for(n, scale=360)
    has_cls = "classified_fraction" in q and q["classified_fraction"].notna().any()
    ncols = 3 if has_cls else 2
    fig, axes = plt.subplots(1, ncols, figsize=(4.2 * ncols + 2, max(3.5, 0.2 * n + 1.8)), sharey=True,
                             gridspec_kw={"width_ratios": [1.3, 1, 0.8][:ncols]})
    y = np.arange(n)
    removed_qc = 100 - q["qc_pass_pairs"] / q["raw_pairs"] * 100
    removed_dec = (q["qc_pass_pairs"] - q["clean_pairs"]) / q["raw_pairs"] * 100
    clean = q["clean_pairs"] / q["raw_pairs"] * 100
    ax = axes[0]
    ax.barh(y, clean, color="#2a7ab9", label="clean")
    ax.barh(y, removed_dec, left=clean, color="#e69f00", label="removed: host / contaminant")
    ax.barh(y, removed_qc, left=clean + removed_dec, color="#999999", label="removed: read QC")
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of raw read pairs")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=1, fontsize=8, frameon=False)
    ax.set_yticks(y)
    ax.set_yticklabels(cols, fontsize=fs)
    style_sample_ticks(cols, low, ax.get_yticklabels())
    ax = axes[1]
    colors = ["#c0392b" if s in low else "#2a7ab9" for s in cols]
    ax.barh(y, q["clean_pairs"], color=colors)
    ax.set_xscale("log")
    if min_pairs > 0:  # --min-clean-pairs 0 flags nothing, and 0 is off a log axis (matplotlib blows the figure up)
        ax.axvline(min_pairs, color="#c0392b", linestyle="--", linewidth=1)
        ax.text(min_pairs, -0.9, f" min {short_count(min_pairs)}", color="#c0392b", fontsize=8, va="bottom")
    for mark in depth_marks:
        ax.axvline(mark, color="#555555", linestyle=":", linewidth=1)
        ax.text(mark, -0.9, f" {short_count(mark)}", color="#555555", fontsize=8, va="bottom")
    ax.set_xlabel("clean read pairs (log)")
    if has_cls:
        ax = axes[2]
        pct = q["classified_fraction"] * 100
        ax.barh(y, pct, color="#56a35a")
        for yi, v in zip(y, pct):
            if pd.isna(v):
                continue
            inside = v > 85
            ax.text(v - 1 if inside else v + 1, yi, f"{v:.0f}%", va="center", ha="right" if inside else "left",
                    fontsize=max(4, fs - 1), color="white" if inside else "#333333")
        ax.set_xlim(0, 100)
        ax.set_xlabel("% classified (profiler)")
    axes[0].invert_yaxis()  # shared y: first sample on top
    for ax in axes:
        for _, start, end in group_spans(cols, groups)[:-1]:
            ax.axhline(end + 0.5, color="#444444", linewidth=0.8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    if any(groups):
        for g, start, end in group_spans(cols, groups):
            axes[-1].text(1.02, (start + end) / 2, str(g), transform=axes[-1].get_yaxis_transform(),
                          va="center", fontsize=fs + 1)
    fig.suptitle("Read fate per sample", fontsize=12)
    fig.tight_layout()
    fw.save(fig, "qc_read_fate", "Read fate per sample",
            "Left: share of raw read pairs kept (blue), removed by host/contaminant filtering (orange) and by "
            "read QC (grey). Middle: clean read pairs on a log scale; red bars and labels mark samples below the "
            "low-depth threshold (red dashed line)"
            + (f"; dotted lines mark {', '.join(short_count(m) for m in depth_marks)} clean read pairs" if depth_marks else "")
            + ". Right: % of reads the taxonomic profiler classified.")


def pick_taxa(wide, per_sample_k, cap=None):
    """Union of each sample's per_sample_k most abundant taxa (unassigned excluded), by maximum abundance."""
    chosen = set()
    for s in wide.columns:
        top = wide[s][wide[s] > 0].nlargest(per_sample_k).index
        chosen.update(top)
    chosen = [t for t in wide.index if t in chosen and not t.startswith("unassigned")]
    ranked = wide.loc[chosen].max(axis=1).sort_values(ascending=False, kind="stable")
    return list(ranked.index[:cap] if cap else ranked.index)


def display_names(names):
    disp = [strip_gtdb(n) for n in names]
    dup = {d for d in disp if disp.count(d) > 1}
    return [n if d in dup else d for n, d in zip(names, disp)]


def parent_of(tx, rank, parent):
    """{taxon name at rank: its parent-rank name} from the lineage table."""
    tax = tx["tax"]
    if parent not in tax.columns or rank not in tax.columns:
        return {}
    pairs = tax[[rank, parent]].replace("", np.nan).dropna().drop_duplicates(rank)
    return dict(zip(pairs[rank], pairs[parent]))


def lineage_order(taxa, tx, rank):
    """Taxa grouped by phylum (largest phylum first), alphabetical within each phylum (GTDB suffixes ignored)."""
    phylum = next((r for r in tx["ranks"] if r.lower() == "phylum"), None)
    to_phylum = parent_of(tx, rank, phylum) if phylum else {}
    ph = {t: to_phylum.get(t, "unassigned") for t in taxa}
    ph_size = pd.Series(ph).value_counts()
    ordered = sorted(taxa, key=lambda t: (-ph_size[ph[t]], ph[t], strip_gtdb(t).lower(), t))
    return ordered, [ph[t] for t in ordered], phylum or "phylum"


def fig_taxonomy(fw, tx, cols, groups, low, per_sample, label):
    sp, ge = species_rank(tx), genus_rank(tx)
    for rank, key in ((sp, "taxonomy_species_heatmap"), (ge, "taxonomy_genus_heatmap")):
        wide = tx["wide"][rank].reindex(columns=cols).fillna(0)
        taxa = pick_taxa(wide, per_sample)
        if not taxa:
            continue
        taxa, phyla, phylum_label = lineage_order(taxa, tx, rank)
        heatmap(fw, key, f"{rank} relative abundance, top {per_sample} per sample ({label})",
                f"Union of each sample's {per_sample} most abundant {rank.lower()} taxa ({len(taxa)} row"
                f"{'s' if len(taxa) != 1 else ''}), so "
                f"every sample's dominant taxa are shown. Rows are grouped by {phylum_label.lower()} (colour strip, "
                f"largest group first) and alphabetical within each. Colour = log10 relative abundance (%); "
                "grey = not detected; "
                "red sample labels mark low-depth samples.",
                wide.loc[taxa], display_names(taxa), cols, groups, low, set(), "log10 relative abundance (%)",
                vmin=-2, vmax=2, row_groups=phyla, row_group_title=phylum_label)
    wide = tx["wide"][ge].reindex(columns=cols).fillna(0)
    taxa = pick_taxa(wide, 3, 12)
    if taxa:
        plot = wide.loc[taxa]
        other = (100 - plot.sum(axis=0)).clip(lower=0)
        n = len(cols)
        fig, ax = plt.subplots(figsize=(min(60, max(8, 0.22 * n + 5)), 5.5))
        colors = category_colors(taxa)
        bottom = np.zeros(n)
        for t, disp in zip(taxa, display_names(taxa)):
            ax.bar(range(n), plot.loc[t].values, bottom=bottom, color=colors[t], width=0.9,
                   label=disp, linewidth=0)
            bottom += plot.loc[t].values
        ax.bar(range(n), other.values, bottom=bottom, color="#d9d9d9", width=0.9, label="Other", linewidth=0)
        ax.set_xticks(range(n))
        ax.set_xticklabels(cols, rotation=90, fontsize=fontsize_for(n))
        style_sample_ticks(cols, low, ax.get_xticklabels())
        ax.set_ylim(0, 100)
        ax.set_ylabel("relative abundance (%)")
        for _, start, end in group_spans(cols, groups)[:-1]:
            ax.axvline(end + 0.5, color="#444444", linewidth=0.8)
        if any(groups):
            label_groups(ax, cols, groups, 101, 5, 8, (min(60, max(8, 0.22 * n + 5)) - 3.5) / n)
        ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        fig.suptitle(f"{ge} composition ({label})", fontsize=12)
        fig.tight_layout()
        fw.save(fig, "taxonomy_genus_bars", f"{ge} composition",
                f"Stacked {ge.lower()} abundance on a linear scale: the {len(taxa)} {ge.lower()} taxa are the union of "
                f"each sample's 3 most abundant (at most 12, by maximum abundance), each with one colour; everything "
                f"else is 'Other'. Unlike the log-scale heatmap, the bars show proportions: which "
                f"{ge.lower()} dominates a sample and how much of the community the named taxa cover.")


def fig_watchlist(fw, wl_wide, wl_long, cols, groups, low, source=""):
    if wl_wide is None or wl_wide.empty:
        return
    det = wl_wide[wl_wide["n_samples_detected"] > 0]
    if det.empty:
        return
    n, m = len(cols), len(det)
    fig, ax = plt.subplots(figsize=(min(60, max(7, 0.2 * n + 4)), max(2.5, 0.3 * m + 2)))
    for i in range(m):
        vals = det.iloc[i][cols].astype(float).values
        idx = np.where(vals > 0)[0]
        if len(idx) == 0:
            continue
        sizes = 12 + 30 * (np.log10(vals[idx]) + 3).clip(0, 5)
        unc = set()
        if wl_long is not None and not wl_long.empty:
            sub = wl_long[(wl_long["taxon"] == det.iloc[i]["taxon"]) & wl_long["uncertain"]]
            unc = set(sub["sample"])
        face = ["none" if cols[j] in unc else "#c0392b" for j in idx]
        ax.scatter(idx, [i] * len(idx), s=sizes, facecolors=face, edgecolors="#c0392b", linewidths=0.8)
    ax.set_yticks(range(m))
    ax.set_yticklabels([f"{t} ({r})" if r == "genus" else t for t, r in zip(det["taxon"], det["rank"])], fontsize=8)
    ax.set_ylim(m - 0.5, -0.5)
    ax.set_xticks(range(n))
    ax.set_xticklabels(cols, rotation=90, fontsize=fontsize_for(n))
    style_sample_ticks(cols, low, ax.get_xticklabels())
    ax.set_xlim(-0.5, n - 0.5)
    for _, start, end in group_spans(cols, groups)[:-1]:
        ax.axvline(end + 0.5, color="#bbbbbb", linewidth=0.8)
    ax.grid(axis="y", color="#eeeeee")
    for v, lab in ((0.01, "0.01%"), (1, "1%"), (50, "50%")):
        ax.scatter([], [], s=12 + 30 * (np.log10(v) + 3), color="#c0392b", label=lab)
    ax.legend(title="abundance", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    ax.set_title("Watchlist taxa", fontsize=12)
    fig.tight_layout()
    fw.save(fig, "watchlist", "Watchlist taxa",
            f"Watchlist taxa found in at least one sample ({len(det)} of the {len(wl_wide)} listed in {source or 'the watchlist'}; "
            "the default list holds the ESKAPE pathogens and other WHO and CDC priority bacteria, see Methods for "
            "each taxon's source). Dot area scales with log abundance. Hollow dots: sylph's ANI interval reaches "
            "below 95%, so the species-level call is uncertain. Detection is not culture: low abundance can be "
            "index hopping or a related strain.")


def fig_args(fw, arg, tables, samples, cols, groups, low, status, top_n, method="K-MARVEL"):
    gene = tables["arg_gene_rpkm"]
    num = gene[cols].apply(pd.to_numeric, errors="coerce")
    hatch = {s for s in cols if status[s] != "ok"}
    sel = gene["tier"] <= 2
    t1 = list(gene.index[(gene["tier"] == 1) & (num.max(axis=1) > 0)])
    t2 = num[sel & (gene["tier"] == 2)]
    t2 = list((t2 > 0).sum(axis=1).sort_values(ascending=False, kind="stable").index[: max(0, top_n - len(t1))])
    rows = t1 + t2
    if rows:
        labels = [f"{g}  [T{int(gene.at[g, 'tier'])}]" for g in rows]
        heatmap(fw, "arg_priority_heatmap", "Acquired ARGs: tier 1 (all) and the most prevalent tier 2 genes",
                "log10 RPKM (supporting read pairs per kb of reference gene per million clean read pairs), one value "
                "per gene (K-MARVEL: its best-supported allele; RGI: summed over the gene's AROs), from "
                f"{method}. Grey = not detected; hatched columns = low depth or {method} not run, where "
                "'not detected' means unknown. Red labels mark low-depth samples.",
                num.loc[rows], labels, cols, groups, low, hatch, "log10 RPKM", cmap="magma_r")
    dc = tables.get("arg_drug_class_rpkm")
    if dc is not None and not dc.empty:
        acq = dc.xs("acquired (tiers 1-2)", level="tier_group") if "acquired (tiers 1-2)" in dc.index.get_level_values(0) else None
        if acq is not None and not acq.empty:
            acqn = acq[cols].apply(pd.to_numeric, errors="coerce")
            acqn = acqn.loc[acqn.max(axis=1).sort_values(ascending=False, kind="stable").index]
            heatmap(fw, "arg_drug_class_heatmap", "Acquired ARG abundance by drug class",
                    "Sum of RPKM of acquired genes (tiers 1-2) conferring resistance to each class. A gene that "
                    "confers resistance to several classes counts in each, so columns are not additive.",
                    acqn, list(acqn.index), cols, groups, low, hatch, "log10 RPKM", cmap="magma_r")
    # load versus depth
    if "clean_pairs" in samples:
        acq_calls = arg[arg["tier"] <= 2]
        n_genes = acq_calls.groupby("sample")["gene"].nunique().reindex(cols).fillna(0)
        load = acq_calls.groupby("sample")["gene_rpkm"].sum().reindex(cols)
        depth = samples.loc[cols, "clean_pairs"].astype(float)
        ran = [s for s in cols if status[s] != NA_NOT_RUN]
        if ran:
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            c = ["#c0392b" if s in low else "#2a7ab9" for s in ran]
            axes[0].scatter(depth[ran], n_genes[ran], c=c, s=22)
            axes[0].set_xscale("log")
            axes[0].set_xlabel("clean read pairs (log)")
            axes[0].set_ylabel("distinct acquired genes (tiers 1-2)")
            axes[1].scatter(depth[ran], load[ran].fillna(0) + 1e-3, c=c, s=22)
            axes[1].set_xscale("log")
            axes[1].set_yscale("log")
            axes[1].set_xlabel("clean read pairs (log)")
            axes[1].set_ylabel("total acquired RPKM (log)")
            for ax in axes:
                for s in ("top", "right"):
                    ax.spines[s].set_visible(False)
            fig.suptitle("ARG detection versus sequencing depth", fontsize=12)
            fig.tight_layout()
            fw.save(fig, "arg_vs_depth", "ARG detection versus sequencing depth",
                    "If the number of genes found rises with depth (left), shallow samples under-report ARGs; "
                    "compare presence/absence only between samples of similar depth. RPKM (right) is depth-"
                    "normalised. Red = below the low-depth threshold.")


TIER_COLORS = {1: "#c0392b", 2: "#e69f00", 3: "#56b4e9", 4: "#999999"}
METHOD_NAMES = {"kmarvel": "K-MARVEL", "rgi": "RGI"}


def fig_method_upset(fw, comp, gene_sets, methods):
    """UpSet plot of the distinct genes each read-level ARG method detected (in any sample), bars split by tier."""
    if not gene_sets:
        return
    tier = comp.drop_duplicates("gene").set_index("gene")["tier"]
    inter = {}
    for g, ms in gene_sets.items():
        inter.setdefault(ms, []).append(g)
    combos = sorted(inter, key=lambda c: (-len(inter[c]), -len(c)))
    names = [METHOD_NAMES.get(m, m) for m in methods]
    nc = len(combos)
    fig = plt.figure(figsize=(max(5.5, 1.1 * nc + 3.5), 2.2 + 0.45 * len(methods) + 2.6), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.3, max(2.0, 0.8 * nc)], height_ratios=[3, 0.45 * len(methods) + 0.2])
    bar = fig.add_subplot(gs[0, 1])
    dots = fig.add_subplot(gs[1, 1], sharex=bar)
    sizes = fig.add_subplot(gs[1, 0], sharey=dots)
    bottom = np.zeros(nc)
    for t in (1, 2, 3, 4):
        h = np.array([sum(1 for g in inter[c] if tier.get(g) == t) for c in combos])
        if h.any():
            bar.bar(range(nc), h, bottom=bottom, color=TIER_COLORS[t], width=0.6, label=f"tier {t}")
            bottom += h
    for x, total in enumerate(bottom):
        bar.text(x, total, f"{int(total)}", ha="center", va="bottom", fontsize=9)
    bar.set_ylabel("distinct genes")
    bar.legend(fontsize=8, frameon=False, loc="upper right")
    bar.tick_params(labelbottom=False, length=0)
    for sp in ("top", "right"):
        bar.spines[sp].set_visible(False)
    for x, c in enumerate(combos):
        ys = [i for i, m in enumerate(methods) if m in c]
        dots.scatter([x] * len(methods), range(len(methods)), s=60, color="#dddddd", zorder=1)
        dots.scatter([x] * len(ys), ys, s=60, color="#333333", zorder=2)
        if len(ys) > 1:
            dots.plot([x, x], [min(ys), max(ys)], color="#333333", linewidth=2, zorder=2)
    dots.set_ylim(len(methods) - 0.5, -0.5)
    dots.set_xlim(-0.6, nc - 0.4)
    dots.set_xticks([])
    dots.tick_params(labelleft=False, length=0)
    for sp in dots.spines.values():
        sp.set_visible(False)
    totals = [sum(1 for ms in gene_sets.values() if m in ms) for m in methods]
    sizes.barh(range(len(methods)), totals, color="#555555", height=0.5)
    for i, v in enumerate(totals):
        sizes.text(v, i, f" {v}", va="center", ha="left", fontsize=8)
    sizes.invert_xaxis()
    sizes.set_yticks(range(len(methods)))
    sizes.set_yticklabels(names, fontsize=9)
    sizes.yaxis.tick_right()
    sizes.set_xlabel("genes per method", fontsize=8)
    for sp in ("top", "left"):
        sizes.spines[sp].set_visible(False)
    title = "Distinct ARGs found by each read-level method"
    fig.suptitle(title, fontsize=12)
    fw.save(fig, "arg_method_upset", title,
            "UpSet plot: each column is one combination of methods (filled dots) and its bar counts the distinct genes "
            "(alleles collapsed, e.g. KPC-2 and KPC-97 are both KPC) that exactly those methods detected in at least "
            "one sample, coloured by tier. Left: genes per method. RGI counts a gene when reads cover at least the "
            "--rgi-min-coverage share of its best reference sequence; only samples every method ran on are compared.")


def fig_method_by_category(fw, comp, methods):
    """Distinct genes per tier/category: found by every method, or by one only (and, for K-MARVEL-only genes,
    whether RGI saw reads below its coverage cutoff)."""
    if comp.empty or len(methods) < 2:
        return
    flags = comp.groupby("gene").agg(tier=("tier", "first"), category=("category", "first"),
                                    **{m: (f"detected_{m}", "any") for m in methods})
    if "rgi_below_cutoff" in comp:
        flags["below"] = comp.groupby("gene")["rgi_below_cutoff"].any()
    else:
        flags["below"] = False
    flags["label"] = "T" + flags["tier"].astype(str) + " " + flags["category"]
    allm = flags[methods].all(axis=1)
    parts = [("all", allm, "#333333", None)]
    colors = ["#2a7ab9", "#e69f00", "#56a35a", "#cc79a7"]
    for i, m in enumerate(methods):
        only = flags[m] & ~flags[[x for x in methods if x != m]].any(axis=1)
        if m == "kmarvel" and "rgi" in methods:
            parts.append((f"{METHOD_NAMES[m]} only, RGI below cutoff", only & flags["below"], colors[i], "///"))
            parts.append((f"{METHOD_NAMES[m]} only, RGI no reads", only & ~flags["below"], colors[i], None))
        else:
            parts.append((f"{METHOD_NAMES.get(m, m)} only", only, colors[i], None))
    parts[0] = (" + ".join(METHOD_NAMES.get(m, m) for m in methods), allm, "#333333", None)
    cats = flags.groupby("label").agg(tier=("tier", "first"), n=("tier", "size"))
    cats = cats.sort_values(["tier", "n"], ascending=[True, False])
    labels = list(cats.index)
    n = len(labels)
    fig, ax = plt.subplots(figsize=(10, max(3, 0.28 * n + 1.8)))
    left = np.zeros(n)
    for name, mask, color, hatch in parts:
        vals = np.array([int((mask & (flags["label"] == lab)).sum()) for lab in labels])
        ax.barh(range(n), vals, left=left, color=color if hatch is None else "white", edgecolor=color,
                hatch=hatch, label=name, height=0.7, linewidth=0.8 if hatch else 0)
        left += vals
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=8)
    for lab, t in zip(ax.get_yticklabels(), cats["tier"]):
        lab.set_color(TIER_COLORS.get(t, "#333333"))
    ax.invert_yaxis()
    ax.set_xlabel("distinct genes")
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    title = "ARG method agreement by tier and category"
    ax.set_title(title, fontsize=12)
    fig.tight_layout()
    fw.save(fig, "arg_method_by_category", title,
            "Distinct genes per reporting tier and category, split by which methods found them in at least one "
            "sample. Hatched: K-MARVEL calls where RGI had reads below the coverage cutoff (partial coverage, not "
            "absence). Category labels are coloured by tier.")


def fig_sarg(fw, mat, cols, groups, low, hatch):
    """ARGs-OAP copies per cell per SARG type, samples as columns."""
    if mat is None or mat.empty:
        return
    num = mat[cols].apply(pd.to_numeric, errors="coerce")
    num = num[num.fillna(0).gt(0).any(axis=1)]
    if num.empty:
        return
    num = num.loc[num.max(axis=1).sort_values(ascending=False, kind="stable").index]
    labels = [str(t).replace("_", " ") for t in num.index]
    heatmap(fw, "arg_sarg_type_heatmap", "ARGs-OAP (SARG): resistance gene copies per cell by type",
            "log10 copies per cell: SARG copies (reads weighted by the share of the reference they cover) summed per "
            "type, divided by the ARGs-OAP cell estimate. SARG types are not CARD drug classes, and copies per cell "
            "is not RPKM: do not compare with the K-MARVEL / RGI figures. Grey = not detected; hatched = low depth, "
            "ARGs-OAP not run or no cell estimate. Red labels mark low-depth samples.",
            num, labels, cols, groups, low, hatch, "log10 copies per cell", cmap="magma_r")


def fig_cog(fw, mat, unit, cols, groups, low, hatch, backend):
    """COG functional categories x samples from the read functional profile."""
    unit_txt = "copies per genome equivalent (CPGE)" if unit == "cpge" else unit
    heatmap(fw, "function_cog_heatmap", f"Functional profile ({backend}): COG categories",
            f"log10 {unit_txt} of each COG functional category, from the read-level profile. An ORF in several "
            "categories counts in each, so the categories do not add up to the annotated reads. Grey = not found; "
            "hatched = low depth, no MicrobeCensus estimate or no functional profile.",
            mat, list(mat.index), cols, groups, low, hatch, f"log10 {unit.upper() if unit == 'cpge' else unit}", cmap="viridis")
