"""Figures of bin/results_summary.py (PNG, PDF and SVG through FigureWriter)."""

import math

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from summary_common import NA_NOT_RUN, natural_key, strip_gtdb  # noqa: E402
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


def label_groups(ax, cols, groups, y0, y1, fontsize):
    """Group names above the columns; alternate two heights so names of narrow groups do not overlap."""
    for i, (g, start, end) in enumerate(group_spans(cols, groups)):
        ax.text((start + end) / 2, y0 if i % 2 == 0 else y1, str(g), ha="center", va="bottom",
                fontsize=fontsize, clip_on=False)


def style_sample_ticks(labels_list, low, axis_labels):
    for lab in axis_labels:
        if lab.get_text() in low:
            lab.set_color("#c0392b")


def heatmap(fw, key, title, caption, mat, row_labels, cols, groups, low, hatch_cols, cbar_label,
            vmin=None, vmax=None, cmap="viridis"):
    data = mat.reindex(columns=cols).astype(float).values
    with np.errstate(divide="ignore", invalid="ignore"):
        logd = np.where(data > 0, np.log10(data), np.nan)
    if np.all(np.isnan(logd)):
        return
    nrow, ncol = logd.shape
    fs_x, fs_y = fontsize_for(ncol), fontsize_for(nrow, scale=420)
    width = min(60, max(7, 0.2 * ncol + 4))
    height = min(60, max(3.5, 0.2 * nrow + 2.8))
    fig, ax = plt.subplots(figsize=(width, height))
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
    if any(groups):
        label_groups(ax, cols, groups, -0.8, -0.8 - max(1.0, nrow * 0.035), fs_x + 1)
    ax.set_xticks(range(ncol))
    ax.set_xticklabels(cols, rotation=90, fontsize=fs_x)
    style_sample_ticks(cols, low, ax.get_xticklabels())
    ax.set_yticks(range(nrow))
    ax.set_yticklabels(row_labels, fontsize=fs_y)
    ax.set_title(title, fontsize=12, pad=38 if any(groups) else 8)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label(cbar_label, fontsize=9)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fw.save(fig, key, title, caption)


def fig_qc(fw, qc, cols, groups, low, min_pairs):
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
    ax.axvline(min_pairs, color="#c0392b", linestyle="--", linewidth=1)
    ax.text(min_pairs, -0.9, f" min {min_pairs:,.0f}", color="#c0392b", fontsize=8, va="bottom")
    ax.set_xlabel("clean read pairs (log)")
    if has_cls:
        ax = axes[2]
        ax.barh(y, q["classified_fraction"] * 100, color="#56a35a")
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
            "low-depth threshold (dashed). Right: fraction of reads the taxonomic profiler classified.")


def pick_taxa(wide, per_sample_k, cap):
    chosen = set()
    for s in wide.columns:
        top = wide[s][wide[s] > 0].nlargest(per_sample_k).index
        chosen.update(top)
    chosen = [t for t in wide.index if t in chosen and not t.startswith("unassigned")]
    ranked = wide.loc[chosen].max(axis=1).sort_values(ascending=False)
    return list(ranked.index[:cap])


def display_names(names):
    disp = [strip_gtdb(n) for n in names]
    dup = {d for d in disp if disp.count(d) > 1}
    return [n if d in dup else d for n, d in zip(names, disp)]


def fig_taxonomy(fw, tx, cols, groups, low, top_n, label):
    sp, ge = species_rank(tx), genus_rank(tx)
    for rank, k, key in ((sp, 5, "taxonomy_species_heatmap"), (ge, 5, "taxonomy_genus_heatmap")):
        wide = tx["wide"][rank].reindex(columns=cols).fillna(0)
        taxa = pick_taxa(wide, k, top_n)
        if not taxa:
            continue
        mat = wide.loc[taxa]
        heatmap(fw, key, f"{rank}: relative abundance ({label})",
                f"Union of each sample's {k} most abundant {rank.lower()} taxa (up to {top_n}), so every sample's "
                f"dominant taxa are shown. Colour = log10 relative abundance (%); grey = not detected. Rows are sorted by "
                f"maximum abundance; red sample labels mark low-depth samples.",
                mat, display_names(taxa), cols, groups, low, set(), "log10 relative abundance (%)",
                vmin=-2, vmax=2)
    wide = tx["wide"][ge].reindex(columns=cols).fillna(0)
    taxa = pick_taxa(wide, 3, 19)
    if taxa:
        plot = wide.loc[taxa]
        other = (100 - plot.sum(axis=0)).clip(lower=0)
        n = len(cols)
        fig, ax = plt.subplots(figsize=(min(60, max(8, 0.22 * n + 5)), 5.5))
        colors = list(plt.get_cmap("tab20").colors)
        bottom = np.zeros(n)
        for i, (t, disp) in enumerate(zip(taxa, display_names(taxa))):
            ax.bar(range(n), plot.loc[t].values, bottom=bottom, color=colors[i % 20], width=0.9,
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
            label_groups(ax, cols, groups, 101, 105, 8)
        ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_title(f"{ge} composition ({label})", fontsize=12, pad=28 if any(groups) else 6)
        fig.tight_layout()
        fw.save(fig, "taxonomy_genus_bars", f"{ge} composition",
                f"Stacked {ge.lower()} abundance; the taxa are the union of each sample's 3 most abundant "
                f"{ge.lower()} taxa (at most 19), so 'Other' stays small. The heatmaps show more taxa.")


def fig_watchlist(fw, wl_wide, wl_long, cols, groups, low):
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
            "Watchlist taxa (priority pathogens by default) per sample; dot area scales with log abundance. "
            "Hollow dots: sylph's ANI interval reaches below 95%, so the species-level call is uncertain. "
            "Detection is not culture: low abundance can be index hopping or a related strain.")


def fig_args(fw, arg, tables, samples, cols, groups, low, status, top_n, group_by, order_by):
    gene = tables["arg_gene_rpkm"]
    num = gene[cols].apply(pd.to_numeric, errors="coerce")
    hatch = {s for s in cols if status[s] != "ok"}
    sel = gene["tier"] <= 2
    t1 = list(gene.index[(gene["tier"] == 1) & (num.max(axis=1) > 0)])
    t2 = num[sel & (gene["tier"] == 2)]
    t2 = list((t2 > 0).sum(axis=1).sort_values(ascending=False).index[: max(0, top_n - len(t1))])
    rows = t1 + t2
    if rows:
        labels = [f"{g}  [T{int(gene.at[g, 'tier'])}]" for g in rows]
        heatmap(fw, "arg_priority_heatmap", "Acquired ARGs: tier 1 (all) and the most prevalent tier 2 genes",
                "log10 RPKM (supporting read pairs per kb of reference gene per million clean read pairs); alleles are "
                "summed per gene. Grey = not detected; hatched columns = low depth or K-MARVEL not run, where "
                "'not detected' means unknown. Red labels mark low-depth samples.",
                num.loc[rows], labels, cols, groups, low, hatch, "log10 RPKM", cmap="magma_r")
    dc = tables.get("arg_drug_class_rpkm")
    if dc is not None and not dc.empty:
        acq = dc.xs("acquired (tiers 1-2)", level="tier_group") if "acquired (tiers 1-2)" in dc.index.get_level_values(0) else None
        if acq is not None and not acq.empty:
            acqn = acq[cols].apply(pd.to_numeric, errors="coerce")
            acqn = acqn.loc[acqn.max(axis=1).sort_values(ascending=False).index]
            heatmap(fw, "arg_drug_class_heatmap", "Acquired ARG abundance by drug class",
                    "Sum of RPKM of acquired genes (tiers 1-2) conferring resistance to each class. A gene that "
                    "confers resistance to several classes counts in each, so columns are not additive.",
                    acqn, list(acqn.index), cols, groups, low, hatch, "log10 RPKM", cmap="magma_r")
    # load versus depth
    if "clean_pairs" in samples:
        acq_calls = arg[arg["tier"] <= 2]
        n_genes = acq_calls.groupby("sample")["gene"].nunique().reindex(cols).fillna(0)
        load = acq_calls.groupby("sample")["rpkm"].sum().reindex(cols)
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
    if group_by and order_by:
        fig_timeline(fw, arg, samples, cols, status, group_by, order_by)


def fig_timeline(fw, arg, samples, cols, status, group_by, order_by):
    t1 = arg[arg["tier"] == 1]
    if t1.empty:
        return
    order_vals = sorted({samples.at[s, order_by] for s in cols}, key=natural_key)
    xi = {v: i for i, v in enumerate(order_vals)}
    rows, cells = [], []
    for g, members in samples.loc[cols].groupby(group_by, sort=False):
        mem = [s for s in cols if s in members.index]
        sub = t1[t1["sample"].isin(mem)]
        for gene in sorted(sub["gene"].unique(), key=str.lower):
            r = len(rows)
            rows.append(f"{g} | {gene}")
            rpkm = sub[sub["gene"] == gene].groupby("sample")["rpkm"].sum()
            for s in mem:
                cells.append((r, xi[samples.at[s, order_by]], rpkm.get(s, np.nan), status[s]))
    if not rows:
        return
    nr, nc = len(rows), len(order_vals)
    fig, ax = plt.subplots(figsize=(max(5, 0.45 * nc + 3.5), max(2.5, 0.28 * nr + 1.6)))
    vals = [math.log10(v) for _, _, v, _ in cells if v and v > 0]
    vmin, vmax = (min(vals), max(vals)) if vals else (0, 1)
    if vmin == vmax:
        vmin, vmax = vmin - 0.5, vmax + 0.5
    cm = plt.get_cmap("magma_r")
    for r, x, v, st in cells:
        if v and v > 0:
            ax.add_patch(Rectangle((x - 0.45, r - 0.4), 0.9, 0.8,
                                   color=cm((math.log10(v) - vmin) / (vmax - vmin) * 0.85 + 0.15)))
        elif st != "ok":
            ax.add_patch(Rectangle((x - 0.45, r - 0.4), 0.9, 0.8, fill=False, hatch="///", edgecolor="#999999"))
        else:
            ax.add_patch(Rectangle((x - 0.45, r - 0.4), 0.9, 0.8, fill=False, edgecolor="#bbbbbb"))
    ax.set_xlim(-0.6, nc - 0.4)
    ax.set_ylim(nr - 0.5, -0.5)
    ax.set_xticks(range(nc))
    ax.set_xticklabels(order_vals, rotation=90 if nc > 12 else 0, fontsize=8)
    ax.set_yticks(range(nr))
    ax.set_yticklabels(rows, fontsize=8)
    ax.set_xlabel(order_by)
    prev = None
    for i, lab in enumerate(rows):
        g = lab.split(" | ")[0]
        if prev is not None and g != prev:
            ax.axhline(i - 0.5, color="#444444", linewidth=0.8)
        prev = g
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(vmin, vmax))
    fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.01).set_label("log10 RPKM", fontsize=9)
    ax.set_title(f"Tier 1 ARGs per {group_by} across {order_by}", fontsize=12)
    fig.tight_layout()
    fw.save(fig, "arg_tier1_timeline", f"Tier 1 ARGs per {group_by} across {order_by}",
            f"Each row is a tier 1 gene found at least once in that {group_by}. Filled = detected (colour = log10 "
            f"RPKM); empty box = sequenced and not detected; hatched = low depth or not run (unknown); no box = no "
            f"sample at that {order_by}.")
