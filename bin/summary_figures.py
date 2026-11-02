"""Figures of bin/results_summary.py, in the shared style of bin/figure_style.py: each figure is written as PDF,
SVG and PNG (600 dpi) at an exact journal width, with its plotted values as CSV. FigureWriter keeps one record per
figure (title, caption, report section, source data and, for the large figures, the data of their interactive
Vega-Lite version, built by bin/summary_interactive.py for the HTML report)."""

import math

import numpy as np
import pandas as pd

import figure_style as fs

fs.setup()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from summary_args import METHOD_LABELS  # noqa: E402
from summary_common import NA_NOT_RUN, strip_gtdb  # noqa: E402
from summary_taxonomy import genus_rank, species_rank  # noqa: E402


class FigureWriter:
    def __init__(self, outdir, min_pairs=None):
        self.dir = outdir / "figures"
        self.dir.mkdir(parents=True, exist_ok=True)
        self.records = []
        self.min_pairs = min_pairs

    def save(self, fig, key, title, caption, source, section, interactive=None):
        """Write the figure files and keep its record. interactive: the data of the HTML's Vega-Lite version
        ({'kind': 'heatmap' | 'read_fate', ...}; None for a figure shown static only)."""
        w, h = fig.get_size_inches()
        fs.save(fig, self.dir, key, source)
        self.records.append({"key": key, "title": title, "caption": caption, "section": section,
                             "width_mm": round(w * 25.4, 1), "height_mm": round(h * 25.4, 1),
                             "interactive": interactive})


short_count = fs.short_count  # also the depth classes of samples_qc (results_summary.py)


def display_names(names):
    """GTDB suffixes removed for display, except where two names would then look alike."""
    disp = [strip_gtdb(n) for n in names]
    dup = {d for d in disp if disp.count(d) > 1}
    return [n if d in dup else d for n, d in zip(names, disp)]


def low_legend(fw, hatch_label=None):
    return fs.low_depth_handle(fw.min_pairs, hatch_label)


def heatmap(fw, key, title, caption, mat, row_labels, cols, groups, low, hatch_cols, cbar_label, section,
            vmin=None, vmax=None, row_groups=None, row_group_title="", feature_name="feature", value_name="value",
            italic=False, hatch_label=None):
    """log10 heatmap of mat (features x samples). Up to 60 samples they are columns, above that rows. Low-depth and
    not-run samples (hatch_cols) are hatched; row_groups (one label per feature, e.g. phylum) adds a colour strip
    with its legend. Colour bar and key below the plot."""
    data = mat.reindex(columns=cols).apply(pd.to_numeric, errors="coerce").astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        logd = np.where(data.values > 0, np.log10(data.values), np.nan)
    if np.all(np.isnan(logd)):
        return
    m, n = logd.shape
    lay = fs.sample_layout(n)
    vertical = lay["vertical"]
    ffs = fs.SMALL_PT if m <= 60 else 5
    fcell = ffs * 0.3528 * 1.05
    shown = [fs.clip_label(x, 36) for x in row_labels]
    longest = max((len(x) for x in shown), default=4) * ffs * 0.3528 * 0.55  # feature label length, mm
    width = fs.width_class(n) if not vertical and longest < 45 else "double"
    slong = max((len(fs.label_text(s, low)) for s in cols), default=4) * (lay["fs"] or 5) * 0.3528 * 0.55
    legend_mm = 14 + (4 * math.ceil(len(set(row_groups)) / 4) if row_groups is not None else 0)
    if vertical:
        h_mm = n * lay["cell_mm"] + longest + legend_mm + 10
    else:
        h_mm = m * fcell + (slong if lay["labels"] else 0) + legend_mm + (6 if any(groups) else 0) + 8
    fig = plt.figure(figsize=fs.figsize(width, h_mm), layout="constrained")
    strip = row_groups is not None
    if vertical:
        gs = fig.add_gridspec(2 if strip else 1, 1, height_ratios=[1.6, max(n * lay["cell_mm"], 20)] if strip else [1])
    else:
        gs = fig.add_gridspec(1, 2 if strip else 1, width_ratios=[1.6, 60] if strip else [1])
    ax = fig.add_subplot(gs[-1, -1])
    z = logd.T if vertical else logd
    vmin = np.nanmin(logd) if vmin is None else vmin
    vmax = np.nanmax(logd) if vmax is None else vmax
    if vmin == vmax:
        vmin, vmax = vmin - 0.5, vmax + 0.5
    im = ax.imshow(np.ma.masked_invalid(z), aspect="auto", cmap=fs.viridis(), vmin=vmin, vmax=vmax,
                   interpolation="nearest")
    fs.hatch_samples(ax, cols, hatch_cols, m, "y" if vertical else "x")
    fs.sample_ticks(ax, cols, low, "y" if vertical else "x", lay["fs"])
    style = "italic" if italic else "normal"
    if vertical:
        ax.set_xticks(range(m), shown, rotation=90, fontsize=ffs, fontstyle=style)
        ax.xaxis.tick_top()
    else:
        ax.set_yticks(range(m), shown, fontsize=ffs, fontstyle=style)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fs.group_brackets(ax, cols, groups, "y" if vertical else "x")
    if strip:
        names = list(dict.fromkeys(row_groups))
        pal = fs.category_colors(names)
        idx = np.array([names.index(g) for g in row_groups])
        sax = fig.add_subplot(gs[0, 0], sharex=ax if vertical else None, sharey=None if vertical else ax)
        sax.imshow(idx[None, :] if vertical else idx[:, None], aspect="auto", interpolation="nearest",
                   cmap=ListedColormap([pal[x] for x in names]), vmin=-0.5, vmax=len(names) - 0.5)
        sax.axis("off")
    heat_mm = fs.WIDTH_MM[width] - (0 if vertical else longest + 8)  # about the width of the heatmap itself
    handles = [Patch(facecolor=fs.NOT_DETECTED, label="not detected")]
    if hatch_cols:
        handles.append(low_legend(fw, hatch_label))
    if strip:
        # one legend: the strip's groups, then the key (two legends side by side overlap on a narrow figure)
        handles = [Patch(color=pal[x], label=d) for x, d in zip(names, display_names(names))] + handles
        fig.legend(handles=handles, title=row_group_title or None, loc="outside lower left", ncol=4,
                   handlelength=1.0, handleheight=0.8, columnspacing=0.8, alignment="left")
    cb = fig.colorbar(im, ax=ax, location="bottom", shrink=1.0 if heat_mm < 60 else 0.45, aspect=25,
                      pad=0.01, anchor=(0, 1))
    cb.set_label(cbar_label, fontsize=fs.SMALL_PT)
    cb.ax.tick_params(labelsize=fs.SMALL_PT, width=0.4, length=1.5)
    cb.outline.set_linewidth(0.4)
    if not strip:
        fig.legend(handles=handles, loc="outside lower right", ncol=len(handles), handlelength=1.2,
                   handleheight=0.8, columnspacing=0.8)
    long = data.copy()
    long.index = list(row_labels)
    long.index.name = feature_name
    src = long.reset_index().melt(id_vars=feature_name, var_name="sample", value_name=value_name)
    gmap = dict(zip(cols, groups))
    src.insert(1, "group", src["sample"].map(gmap))
    src["low_depth"] = src["sample"].isin(low)
    src["hatched"] = src["sample"].isin(hatch_cols)
    if strip:
        disp = dict(zip(dict.fromkeys(row_groups), display_names(list(dict.fromkeys(row_groups)))))
        src["row_group"] = src[feature_name].map(dict(zip(row_labels, [disp[g] for g in row_groups])))
    src = src[[feature_name, "sample", "group", "low_depth", "hatched"] + (["row_group"] if strip else [])
              + [value_name]]
    fw.save(fig, key, title, caption, src, section,
            {"kind": "heatmap", "frame": src, "feature": feature_name, "value": value_name,
             "features": list(row_labels), "samples": list(cols), "groups": list(groups), "low": sorted(low),
             "domain": [float(vmin), float(vmax)], "legend": cbar_label, "row_group": "row_group" if strip else None,
             "row_group_title": row_group_title, "italic": italic})


def fig_qc(fw, qc, cols, groups, low, min_pairs, depth_marks=()):
    """Read fate per sample: share of raw read pairs kept / removed, clean read pairs (log), % classified."""
    if "raw_pairs" not in qc or qc["raw_pairs"].isna().all():
        return
    q = qc.reindex(cols)
    n = len(cols)
    lay = fs.sample_layout(n)
    has_cls = "classified_fraction" in q and q["classified_fraction"].notna().any()
    clean = (q["clean_pairs"] / q["raw_pairs"] * 100).values
    host = ((q["qc_pass_pairs"] - q["clean_pairs"]) / q["raw_pairs"] * 100).values
    qcr = (100 - q["qc_pass_pairs"] / q["raw_pairs"] * 100).values
    width = "double" if has_cls or n > 14 else "single"
    h_mm = n * lay["cell_mm"] + 26 + (4 if depth_marks else 0)
    ratios = [1.4, 1, 0.8][: 3 if has_cls else 2]
    fig, axes = plt.subplots(1, len(ratios), figsize=fs.figsize(width, h_mm), sharey=True, layout="constrained",
                             gridspec_kw={"width_ratios": ratios})
    y = np.arange(n)
    a1, a2 = axes[0], axes[1]
    parts = [(clean, "clean", fs.CATEGORY[0]), (host, "removed: host / contaminant", fs.CATEGORY[1]),
             (qcr, "removed: read QC", "#9A9A9A")]
    left = np.zeros(n)
    for vals, lab, col in parts:
        v = np.nan_to_num(vals)
        a1.barh(y, v, left=left, color=col, label=lab, height=0.8, lw=0)
        left += v
    a1.set_xlim(0, 100)
    a1.set_xlabel("% of raw read pairs")
    fs.sample_ticks(a1, cols, low, "y", lay["fs"])
    a1.set_ylim(n - 0.5, -0.5)
    pairs = q["clean_pairs"].astype(float).values
    bars = a2.barh(y, pairs, color=fs.CATEGORY[0], height=0.8, lw=0)
    for b, s in zip(bars, cols):
        if s in low:
            b.set_facecolor("#9A9A9A")
            b.set_hatch("////")
            b.set_edgecolor("white")
    a2.set_xscale("log")
    fs.count_axis(a2.xaxis)
    finite = pairs[np.isfinite(pairs) & (pairs > 0)]
    if len(finite):
        lo = 10 ** math.floor(math.log10(min(finite.min(), min_pairs or finite.min()) / 1.5))
        a2.set_xlim(lo, 10 ** math.ceil(math.log10(max(finite.max(), max(depth_marks, default=0)) * 1.05)))
    if min_pairs > 0:  # --min-clean-pairs 0 flags nothing, and 0 is off a log axis
        a2.axvline(min_pairs, color=fs.LOW_DEPTH, ls="--", lw=0.6)
    for mark in depth_marks:
        a2.axvline(mark, color="#555555", ls=":", lw=0.6)
        a2.text(mark, 1.0, f" {fs.short_count(mark)}", transform=a2.get_xaxis_transform(), color="#555555",
                fontsize=5, va="bottom", ha="left", clip_on=False)
    a2.set_xlabel("clean read pairs (log scale)")
    if has_cls:
        a3 = axes[2]
        pct = (q["classified_fraction"] * 100).values
        a3.barh(y, np.nan_to_num(pct), color=fs.CATEGORY[2], height=0.8, lw=0)
        a3.set_xlim(0, 100)
        a3.set_xlabel("% classified (profiler)")
    for ax in axes:
        ax.tick_params(axis="y", length=0)
        for g, a, b in fs.group_spans(cols, groups)[1:] if any(groups) else []:
            ax.axhline(a - 0.5, color="#333333", lw=0.4)
    fs.group_brackets(axes[-1], cols, groups, "y", separator=None, offset=0.03)
    handles, labels = a1.get_legend_handles_labels()
    if low:
        h = low_legend(fw)
        handles.append(h)
        labels.append(h.get_label() + " (dashed line)")
    fig.legend(handles, labels, loc="outside lower center", ncol=4 if width == "double" else 2, handlelength=1.2,
               columnspacing=1.0)
    src = pd.DataFrame({"sample": cols, "group": groups, "low_depth": [s in low for s in cols],
                        "raw_pairs": q["raw_pairs"].values, "qc_pass_pairs": q["qc_pass_pairs"].values,
                        "clean_pairs": q["clean_pairs"].values, "pct_clean": clean, "pct_removed_host": host,
                        "pct_removed_qc": qcr})
    if has_cls:
        src["pct_classified"] = (q["classified_fraction"] * 100).values
    marks = f"; dotted lines mark {', '.join(fs.short_count(m) for m in depth_marks)} clean read pairs" \
        if depth_marks else ""
    fw.save(fig, "qc_read_fate", "Read fate per sample",
            "Left: share of raw read pairs kept (clean), removed as host or contaminant, and removed by read QC. "
            f"Middle: clean read pairs (log scale); the dashed line is the low-depth cutoff ({min_pairs:,} pairs), "
            f"and hatched grey bars and {fs.DAGGER} labels mark samples below it{marks}."
            + (" Right: % of reads the taxonomic profiler classified." if has_cls else ""),
            src, "overview", {"kind": "read_fate", "frame": src, "samples": list(cols), "groups": list(groups),
                        "low": sorted(low), "min_pairs": min_pairs, "classified": bool(has_cls)})


def pick_taxa(wide, per_sample_k, cap=None):
    """Union of each sample's per_sample_k most abundant taxa (unassigned excluded), by maximum abundance."""
    chosen = set()
    for s in wide.columns:
        top = wide[s][wide[s] > 0].nlargest(per_sample_k).index
        chosen.update(top)
    chosen = [t for t in wide.index if t in chosen and not t.startswith("unassigned")]
    ranked = wide.loc[chosen].max(axis=1).sort_values(ascending=False, kind="stable")
    return list(ranked.index[:cap] if cap else ranked.index)


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
        heatmap(fw, key, f"{rank} relative abundance ({label})",
                f"Union of each sample's {per_sample} most abundant {rank.lower()} taxa ({len(taxa)} row"
                f"{'s' if len(taxa) != 1 else ''}), so every sample's dominant taxa are shown, grouped by "
                f"{phylum_label.lower()} (colour strip, largest group first) and alphabetical within each. Colour: "
                f"log10 relative abundance (%); grey: not detected; {fs.DAGGER} and vermillion labels: low-depth "
                "samples.",
                wide.loc[taxa], display_names(taxa), cols, groups, low, set(), "log10 relative abundance (%)",
                "taxonomy", vmin=-2, vmax=2, row_groups=phyla, row_group_title=phylum_label,
                feature_name=rank.lower(), value_name="relative_abundance_pct", italic=True)
    wide = tx["wide"][ge].reindex(columns=cols).fillna(0)
    taxa = pick_taxa(wide, 3, 12)
    if not taxa:
        return
    plot = wide.loc[taxa]
    other = (100 - plot.sum(axis=0)).clip(lower=0)
    n = len(cols)
    lay = fs.sample_layout(n)
    width = fs.width_class(n)
    fig, ax = plt.subplots(figsize=fs.figsize(width, n * lay["cell_mm"] + 22 + 3 * math.ceil(len(taxa) / 3)),
                           layout="constrained")
    names = display_names(taxa)
    pal = fs.category_colors(names + ["Other"])
    y = np.arange(n)
    left = np.zeros(n)
    for t, disp in zip(taxa, names):
        ax.barh(y, plot.loc[t].values, left=left, color=pal[disp], height=0.85, label=disp, lw=0)
        left += plot.loc[t].values
    ax.barh(y, other.values, left=left, color=pal["Other"], height=0.85, label="Other", lw=0)
    ax.set_xlim(0, 100)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlabel("relative abundance (%)")
    fs.sample_ticks(ax, cols, low, "y", lay["fs"])
    ax.tick_params(axis="y", length=0)
    for g, a, b in fs.group_spans(cols, groups)[1:] if any(groups) else []:
        ax.axhline(a - 0.5, color="white", lw=0.8)
    fs.group_brackets(ax, cols, groups, "y", separator=None, offset=0.02)
    leg = fig.legend(loc="outside lower center", ncol=3, handlelength=1.0, handleheight=0.8, columnspacing=0.8)
    for t in leg.get_texts():
        if t.get_text() != "Other":
            t.set_fontstyle("italic")
    src = plot.copy()
    src.index = names
    src.loc["Other"] = other
    src = src.rename_axis(ge.lower()).reset_index().melt(id_vars=ge.lower(), var_name="sample",
                                                         value_name="relative_abundance_pct")
    src.insert(2, "group", src["sample"].map(dict(zip(cols, groups))))
    fw.save(fig, "taxonomy_genus_bars", f"{ge} composition ({label})",
            f"Stacked {ge.lower()} abundance on a linear scale: the {len(taxa)} {ge.lower()} taxa are the union of "
            f"each sample's 3 most abundant (at most 12, by maximum abundance), each with one colour; everything "
            f"else is 'Other'. Unlike the log-scale heatmap, the bars show proportions: which {ge.lower()} dominates a "
            "sample and how much of the community the named taxa cover.", src, "taxonomy")


def fig_watchlist(fw, wl_wide, wl_long, cols, groups, low, source=""):
    if wl_wide is None or wl_wide.empty:
        return
    det = wl_wide[wl_wide["n_samples_detected"] > 0]
    if det.empty:
        return
    n, m = len(cols), len(det)
    lay = fs.sample_layout(n)
    vertical = lay["vertical"]
    labels = [f"{t} (genus)" if r == "genus" else t for t, r in zip(det["taxon"], det["rank"])]
    width = fs.width_class(n) if not vertical else "double"
    h_mm = (n * lay["cell_mm"] + 45) if vertical else (m * 3.0 + 30 + (6 if any(groups) else 0))
    fig, ax = plt.subplots(figsize=fs.figsize(width, h_mm), layout="constrained")
    rows = []
    for i in range(m):
        vals = det.iloc[i][cols].astype(float).values
        unc = set()
        if wl_long is not None and not wl_long.empty:
            sub = wl_long[(wl_long["taxon"] == det.iloc[i]["taxon"]) & wl_long["uncertain"]]
            unc = set(sub["sample"])
        for j, v in enumerate(vals):
            rows.append({"taxon": labels[i], "sample": cols[j], "group": groups[j], "low_depth": cols[j] in low,
                         "relative_abundance_pct": v, "species_call_uncertain": cols[j] in unc})
        idx = np.where(vals > 0)[0]
        if not len(idx):
            continue
        sizes = 2 + 8 * (np.log10(vals[idx]) + 3).clip(0, 5)
        face = ["none" if cols[j] in unc else fs.LOW_DEPTH for j in idx]
        xs, ys = (np.full(len(idx), i), idx) if vertical else (idx, np.full(len(idx), i))
        ax.scatter(xs, ys, s=sizes, facecolors=face, edgecolors=fs.LOW_DEPTH, linewidths=0.5)
    feat_axis, samp_axis = ("x", "y") if vertical else ("y", "x")
    if vertical:
        ax.set_xticks(range(m), labels, rotation=90, fontsize=fs.SMALL_PT, fontstyle="italic")
        ax.xaxis.tick_top()
        ax.set_xlim(-0.5, m - 0.5)
        ax.set_ylim(n - 0.5, -0.5)
    else:
        ax.set_yticks(range(m), labels, fontsize=fs.SMALL_PT, fontstyle="italic")
        ax.set_ylim(m - 0.5, -0.5)
        ax.set_xlim(-0.5, n - 0.5)
    fs.sample_ticks(ax, cols, low, samp_axis, lay["fs"])
    ax.grid(axis=feat_axis, color="#EEEEEE", lw=0.4)
    ax.tick_params(length=0)
    fs.group_brackets(ax, cols, groups, samp_axis, separator="#CCCCCC")
    handles = [Line2D([], [], ls="", marker="o", markersize=math.sqrt(2 + 8 * (math.log10(v) + 3)),
                      color=fs.LOW_DEPTH, label=lab) for v, lab in ((0.01, "0.01 %"), (1, "1 %"), (50, "50 %"))]
    handles.append(Line2D([], [], ls="", marker="o", markersize=3, markerfacecolor="none", color=fs.LOW_DEPTH,
                          label="species call uncertain"))
    fig.legend(handles=handles, title="relative abundance", loc="outside lower center", ncol=4, columnspacing=1.0)
    fw.save(fig, "watchlist", "Watchlist taxa",
            f"Watchlist taxa found in at least one sample ({len(det)} of the {len(wl_wide)} listed in "
            f"{source or 'the watchlist'}; the default list holds the ESKAPE pathogens and other WHO and CDC priority "
            "bacteria, see Methods for each taxon's source). Dot area scales with log abundance. Hollow dots: sylph's "
            "ANI interval reaches below 95 %, so the species-level call is uncertain. Detection is not culture: low "
            "abundance can be index hopping or a related strain.", pd.DataFrame(rows), "taxonomy")


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
        labels = [f"{g} [T{int(gene.at[g, 'tier'])}]" for g in rows]
        heatmap(fw, "arg_priority_heatmap", "Acquired ARGs: tier 1 and the most prevalent tier 2 genes",
                "log10 RPKM (supporting read pairs per kb of reference gene per million clean read pairs), one value "
                "per gene (K-MARVEL: its best-supported allele; RGI: summed over the gene's AROs), from "
                f"{method}; every tier 1 gene found, then the tier 2 genes found in the most samples (at most "
                f"{top_n} rows). Grey: not detected; hatched: low depth or {method} not run, where 'not detected' "
                "means unknown.",
                num.loc[rows], labels, cols, groups, low, hatch, "log10 RPKM", "args",
                feature_name="gene", value_name="rpkm", hatch_label=f"{fs.DAGGER} low depth or not run")
    dc = tables.get("arg_drug_class_rpkm")
    if dc is not None and not dc.empty and "acquired (tiers 1-2)" in dc.index.get_level_values(0):
        acq = dc.xs("acquired (tiers 1-2)", level="tier_group")
        if not acq.empty:
            acqn = acq[cols].apply(pd.to_numeric, errors="coerce")
            acqn = acqn.loc[acqn.max(axis=1).sort_values(ascending=False, kind="stable").index]
            heatmap(fw, "arg_drug_class_heatmap", "Acquired ARG abundance by drug class",
                    "Sum of the RPKM of the acquired genes (tiers 1-2) conferring resistance to each class. A gene "
                    "that confers resistance to several classes counts in each, so the rows are not additive. Grey: "
                    "not detected; hatched: low depth or not run.",
                    acqn, list(acqn.index), cols, groups, low, hatch, "log10 RPKM", "args",
                    feature_name="drug_class", value_name="rpkm", hatch_label=f"{fs.DAGGER} low depth or not run")
    if "clean_pairs" not in samples:
        return
    acq_calls = arg[arg["tier"] <= 2]
    n_genes = acq_calls.groupby("sample")["gene"].nunique().reindex(cols).fillna(0)
    load = acq_calls.groupby("sample")["gene_rpkm"].sum().reindex(cols)
    depth = samples.loc[cols, "clean_pairs"].astype(float)
    ran = [s for s in cols if status[s] != NA_NOT_RUN]
    if not ran:
        return
    fig, axes = plt.subplots(1, 2, figsize=fs.figsize("single", 48), layout="constrained")
    c = [fs.LOW_DEPTH if s in low else fs.CATEGORY[0] for s in ran]
    axes[0].scatter(depth[ran], n_genes[ran], c=c, s=6, lw=0)
    axes[0].set_ylabel("distinct acquired genes")
    axes[1].scatter(depth[ran], load[ran].fillna(0) + 1e-3, c=c, s=6, lw=0)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("total acquired RPKM (log)")
    for ax in axes:
        ax.set_xscale("log")
        fs.count_axis(ax.xaxis)
        ax.set_xlabel("clean read pairs (log)")
    if low:
        fig.legend(handles=[Line2D([], [], ls="", marker="o", markersize=2.5, color=fs.CATEGORY[0], label="sample"),
                            Line2D([], [], ls="", marker="o", markersize=2.5, color=fs.LOW_DEPTH,
                                   label=f"low depth ({fs.DAGGER})")], loc="outside lower center", ncol=2)
    src = pd.DataFrame({"sample": ran, "group": [groups[cols.index(s)] for s in ran],
                        "low_depth": [s in low for s in ran], "clean_pairs": depth[ran].values,
                        "acquired_genes": n_genes[ran].values, "acquired_rpkm": load[ran].fillna(0).values})
    fw.save(fig, "arg_vs_depth", "ARG detection versus sequencing depth",
            "Acquired genes (tiers 1-2) per sample against its clean read pairs. If the number of genes found rises "
            "with depth (left), shallow samples under-report ARGs: compare presence and absence only between samples "
            "of similar depth. Total RPKM (right) is normalised to depth. Vermillion: below the low-depth threshold.",
            src, "args")


TIER_COLORS = {1: "#D55E00", 2: "#E69F00", 3: "#56B4E9", 4: "#9A9A9A"}


def fig_method_upset(fw, comp, gene_sets, methods):
    """UpSet plot of the distinct genes each read-level ARG method detected (in any sample), bars split by tier."""
    if not gene_sets:
        return
    tier = comp.drop_duplicates("gene").set_index("gene")["tier"]
    inter = {}
    for g, ms in gene_sets.items():
        inter.setdefault(ms, []).append(g)
    combos = sorted(inter, key=lambda c: (-len(inter[c]), -len(c)))
    names = [METHOD_LABELS.get(m, m) for m in methods]
    nc = len(combos)
    fig = plt.figure(figsize=fs.figsize("single", 55 + 4 * len(methods)), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[1, max(1.6, 0.7 * nc)], height_ratios=[3, 0.3 * len(methods) + 0.2])
    bar = fig.add_subplot(gs[0, 1])
    dots = fig.add_subplot(gs[1, 1], sharex=bar)
    sizes = fig.add_subplot(gs[1, 0], sharey=dots)
    bottom = np.zeros(nc)
    rows = []
    for t in (1, 2, 3, 4):
        h = np.array([sum(1 for g in inter[c] if tier.get(g) == t) for c in combos])
        for c, v in zip(combos, h):
            rows.append({"methods": " + ".join(METHOD_LABELS.get(m, m) for m in methods if m in c), "tier": t,
                         "genes": int(v)})
        if h.any():
            bar.bar(range(nc), h, bottom=bottom, color=TIER_COLORS[t], width=0.6, label=f"tier {t}", lw=0)
            bottom += h
    for x, total in enumerate(bottom):
        bar.text(x, total, f"{int(total)}", ha="center", va="bottom", fontsize=fs.SMALL_PT)
    bar.set_ylabel("distinct genes")
    bar.legend(loc="upper right", handlelength=1.0)
    bar.tick_params(labelbottom=False, bottom=False)
    for x, c in enumerate(combos):
        ys = [i for i, m in enumerate(methods) if m in c]
        dots.scatter([x] * len(methods), range(len(methods)), s=14, color="#DDDDDD", zorder=1, lw=0)
        dots.scatter([x] * len(ys), ys, s=14, color="#333333", zorder=2, lw=0)
        if len(ys) > 1:
            dots.plot([x, x], [min(ys), max(ys)], color="#333333", lw=1, zorder=2)
    dots.set_ylim(len(methods) - 0.5, -0.5)
    dots.set_xlim(-0.6, nc - 0.4)
    dots.set_xticks([])
    dots.tick_params(labelleft=False, length=0)
    for sp in dots.spines.values():
        sp.set_visible(False)
    totals = [sum(1 for ms in gene_sets.values() if m in ms) for m in methods]
    sizes.barh(range(len(methods)), totals, color="#555555", height=0.5, lw=0)
    for i, v in enumerate(totals):
        sizes.text(v, i, f"{v} ", va="center", ha="right", fontsize=fs.SMALL_PT, color="white")
    sizes.invert_xaxis()
    sizes.set_yticks(range(len(methods)), names)
    sizes.yaxis.tick_right()
    sizes.tick_params(axis="y", length=0)
    sizes.set_xlabel("genes per method")
    sizes.spines["left"].set_visible(False)
    sizes.spines["right"].set_visible(True)
    fw.save(fig, "arg_method_upset", "Distinct ARGs found by each read-level method",
            "UpSet plot: each column is one combination of methods (filled dots), and its bar counts the distinct genes "
            "(alleles collapsed: KPC-2 and KPC-97 are both KPC) that exactly those methods detected in at least one "
            "sample, coloured by tier. Left: genes per method. RGI counts a gene when reads cover at least the "
            "--rgi-min-coverage share of its best reference sequence; only samples every method ran on are compared.",
            pd.DataFrame(rows), "compare")


def fig_method_by_category(fw, comp, methods):
    """Distinct genes per tier/category: found by every method, or by one only (and, for K-MARVEL-only genes,
    whether RGI saw reads below its coverage cutoff)."""
    if comp.empty or len(methods) < 2:
        return
    flags = comp.groupby("gene").agg(tier=("tier", "first"), category=("category", "first"),
                                    **{m: (f"detected_{m}", "any") for m in methods})
    flags["below"] = comp.groupby("gene")["rgi_below_cutoff"].any() if "rgi_below_cutoff" in comp else False
    flags["label"] = "T" + flags["tier"].astype(str) + " " + flags["category"]
    allm = flags[methods].all(axis=1)
    parts = [(" + ".join(METHOD_LABELS.get(m, m) for m in methods), allm, "#333333", None)]
    colors = [fs.CATEGORY[0], fs.CATEGORY[1], fs.CATEGORY[2], fs.CATEGORY[3]]
    for i, m in enumerate(methods):
        only = flags[m] & ~flags[[x for x in methods if x != m]].any(axis=1)
        if m == "kmarvel" and "rgi" in methods:
            parts.append((f"{METHOD_LABELS[m]} only, RGI below cutoff", only & flags["below"], colors[i], "////"))
            parts.append((f"{METHOD_LABELS[m]} only, RGI no reads", only & ~flags["below"], colors[i], None))
        else:
            parts.append((f"{METHOD_LABELS.get(m, m)} only", only, colors[i], None))
    cats = flags.groupby("label").agg(tier=("tier", "first"), n=("tier", "size"))
    cats = cats.sort_values(["tier", "n"], ascending=[True, False])
    labels = list(cats.index)
    n = len(labels)
    fig, ax = plt.subplots(figsize=fs.figsize("double", n * 2.6 + 24), layout="constrained")
    left = np.zeros(n)
    rows = []
    for name, mask, color, hatch in parts:
        vals = np.array([int((mask & (flags["label"] == lab)).sum()) for lab in labels])
        rows += [{"category": lab, "found_by": name, "genes": int(v)} for lab, v in zip(labels, vals)]
        ax.barh(range(n), vals, left=left, color=color if hatch is None else "white", edgecolor=color,
                hatch=hatch, label=name, height=0.7, lw=0.5 if hatch else 0)
        left += vals
    ax.set_yticks(range(n), labels)
    for lab, t in zip(ax.get_yticklabels(), cats["tier"]):
        lab.set_color(TIER_COLORS.get(t, "#333333") if t != 4 else "#555555")
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlabel("distinct genes")
    fig.legend(loc="outside lower center", ncol=len(parts), handlelength=1.2, columnspacing=1.0)
    fw.save(fig, "arg_method_by_category", "ARG method agreement by tier and category",
            "Distinct genes per reporting tier and category, split by which methods found them in at least one "
            "sample. Hatched: K-MARVEL calls where RGI had reads below the coverage cutoff (partial coverage, not "
            "absence). Category labels are coloured by tier.", pd.DataFrame(rows), "compare")


def fig_sarg(fw, mat, cols, groups, low, hatch):
    """ARGs-OAP copies per cell per SARG type, samples as columns."""
    if mat is None or mat.empty:
        return
    num = mat[cols].apply(pd.to_numeric, errors="coerce")
    num = num[num.fillna(0).gt(0).any(axis=1)]
    if num.empty:
        return
    num = num.loc[num.max(axis=1).sort_values(ascending=False, kind="stable").index]
    heatmap(fw, "arg_sarg_type_heatmap", "ARGs-OAP (SARG): resistance gene copies per cell by type",
            "log10 copies per cell: SARG copies (reads weighted by the share of the reference they cover) summed per "
            "type, divided by the ARGs-OAP cell estimate; types ordered by their highest value. SARG types are not "
            "CARD drug classes, and copies per cell is not RPKM: do not compare with the K-MARVEL / RGI figures. Grey: "
            "not detected; hatched: low depth, ARGs-OAP not run or no cell estimate.",
            num, [str(t).replace("_", " ") for t in num.index], cols, groups, low, hatch, "log10 copies per cell",
            "sarg", feature_name="sarg_type", value_name="copies_per_cell",
            hatch_label=f"{fs.DAGGER} low depth, not run or not normalised")


def fig_cog(fw, mat, unit, cols, groups, low, hatch, backend, key="function_cog_heatmap", section="functional",
            level="read"):
    """COG functional categories x samples (read- or contig-level profile)."""
    num = mat.reindex(columns=cols).apply(pd.to_numeric, errors="coerce")
    if num.fillna(0).le(0).all().all():
        return
    unit_txt = "copies per genome equivalent (CPGE)" if unit == "cpge" else unit.upper() if unit == "tpm" else unit
    desc = mat["description"] if "description" in mat.columns else pd.Series("", index=mat.index)
    labels = [f"{a} {d}".strip() if isinstance(d, str) else str(a) for a, d in zip(mat.index, desc)]
    what = (f"the read-level profile ({backend})" if level == "read" else
            "the contig-level profile (eggNOG-mapper on the predicted genes, featureCounts per sample)")
    heatmap(fw, key, f"Functional profile ({backend}): COG categories" if level == "read"
            else "Contig functional profile: COG categories",
            f"log10 {unit_txt} of each COG functional category, from {what}. A gene in several categories counts in "
            "each, so the categories do not add up. Grey: not found; hatched: low depth, "
            + ("no MicrobeCensus estimate or no functional profile." if level == "read" else
               "no MicrobeCensus estimate or no contig profile."),
            num, labels, cols, groups, low, hatch, f"log10 {'CPGE' if unit == 'cpge' else unit_txt}", section,
            feature_name="cog_category", value_name=unit if unit in ("cpge", "tpm") else "abundance",
            hatch_label=f"{fs.DAGGER} low depth or no estimate")


MIMAG_COLORS = {"high": "#009E73", "medium": "#E69F00", "low": "#9A9A9A"}


def fig_mag_quality(fw, mags):
    """CheckM2 completeness against contamination of every refined MAG, coloured by MIMAG tier."""
    if mags is None or mags.empty:
        return
    m = mags.copy()
    m["completeness_pct"] = pd.to_numeric(m["completeness_pct"], errors="coerce")
    m["contamination_pct"] = pd.to_numeric(m["contamination_pct"], errors="coerce")
    fig, ax = plt.subplots(figsize=fs.figsize("single", 62), layout="constrained")
    for tier in ("low", "medium", "high"):
        sub = m[m["quality_tier"] == tier]
        if len(sub):
            ax.scatter(sub["contamination_pct"], sub["completeness_pct"], s=7, color=MIMAG_COLORS[tier], lw=0,
                       alpha=0.9, label=f"{tier} ({len(sub)})")
    for x, y in ((5, 90), (10, 50)):
        ax.axhline(y, color="#BBBBBB", lw=0.4, ls=":", zorder=0)
        ax.axvline(x, color="#BBBBBB", lw=0.4, ls=":", zorder=0)
    ax.set_xlabel("contamination (%, CheckM2)")
    ax.set_ylabel("completeness (%, CheckM2)")
    ax.set_xlim(left=0)
    ax.set_ylim(top=101)
    fig.legend(title="MIMAG tier", loc="outside right upper")
    src = m[["assembly", "batch", "bin", "completeness_pct", "contamination_pct", "quality_tier", "tier_basis"]]
    n_ass = m["assembly"].nunique()
    fw.save(fig, "mag_quality", "MAG quality",
            f"CheckM2 completeness against contamination of the {len(m)} refined MAGs of {n_ass} assembl"
            f"{'y' if n_ass == 1 else 'ies'}, coloured by MIMAG tier (high: > 90 % complete, < 5 % contaminated, with "
            "the 23S, 16S and 5S rRNA genes and at least 18 tRNA types from Bakta; medium: >= 50 %, < 10 %). Dotted "
            "lines: the MIMAG cutoffs. Without Bakta no MAG is high (tier_basis checkm2_only).", src, "mags")


def fig_mag_abundance(fw, ab, mags, cols, groups, low, max_rows=80):
    """Share of each sample's sequenced bases on each MAG (rows: MAGs, grouped by phylum)."""
    if ab is None or ab.empty or mags is None or mags.empty:
        return
    a = ab.copy()
    a["pct_sample_bases"] = pd.to_numeric(a["pct_sample_bases"], errors="coerce")
    # a MAG is batch + assembly + bin: two merged co-assemblies are both called 'coassembly'
    multi = mags["batch"].nunique() > 1

    def mag_key(d):
        key = d["assembly"].astype(str) + " " + d["bin"].astype(str)
        return d["batch"].astype(str) + ": " + key if multi else key
    a["mag"] = mag_key(a)
    mat = a.pivot_table(index="mag", columns="sample", values="pct_sample_bases", aggfunc="first")
    mat = mat.reindex(columns=[c for c in cols if c in mat.columns])
    if mat.empty or mat.fillna(0).le(0).all().all():
        return
    info = mags.assign(mag=mag_key(mags)).drop_duplicates("mag").set_index("mag")
    top = mat.max(axis=1).sort_values(ascending=False, kind="stable").index[:max_rows]
    ph = {k: (strip_gtdb(str(info.at[k, "phylum"])) if k in info.index and str(info.at[k, "phylum"]) not in
              ("", "nan") else "unclassified") for k in top}
    size = pd.Series(ph).value_counts()
    rows = sorted(top, key=lambda k: (-size[ph[k]], ph[k], -float(mat.loc[k].max())))
    genus = {k: strip_gtdb(str(info.at[k, "genus"])) if k in info.index and str(info.at[k, "genus"]) not in ("", "nan")
             else "" for k in rows}
    labels = [f"{k} {genus[k]}".strip() for k in rows]
    scols = list(mat.columns)
    sgroups = [groups[cols.index(s)] for s in scols]
    capped = len(mat) > max_rows
    heatmap(fw, "mag_abundance_heatmap", "MAG abundance per sample",
            "Share of each sample's sequenced bases on each MAG (pct_sample_bases, log10 colour scale), "
            + (f"the {max_rows} MAGs with the highest share of the {len(mat)}, " if capped else f"{len(mat)} MAGs, ")
            + "grouped by GTDB phylum (colour strip) and labelled by assembly, bin and genus. Grey: no reads of the "
            "sample on the MAG, or (per-sample assembly) a MAG of another sample's assembly.",
            mat.loc[rows], labels, scols, sgroups, low, set(), "log10 % of sample bases", "mags",
            row_groups=[ph[k] for k in rows], row_group_title="phylum", feature_name="mag",
            value_name="pct_sample_bases")


def fig_arg_contig_hosts(fw, hosts, top=20, top_classes=8):
    """ARG genes on contigs per host genus (all assemblies), stacked by ARG class."""
    if hosts is None or hosts.empty:
        return
    h = hosts.copy()
    h["arg_genes"] = pd.to_numeric(h["arg_genes"], errors="coerce").fillna(0)
    h["host"] = np.where(h["host_taxon"].fillna("").astype(str) == "", "unclassified", h["host_taxon"].astype(str))
    h.loc[(h["host_rank"] != "genus") & (h["host"] != "unclassified"), "host"] = \
        h["host"] + " (" + h["host_rank"].astype(str) + ")"
    by_host = h.groupby("host")["arg_genes"].sum().sort_values(ascending=False, kind="stable")
    hosts_kept = [x for x in by_host.index if x != "unclassified"][:top]
    if "unclassified" in by_host.index:
        hosts_kept.append("unclassified")
    classes = h.groupby("arg_class")["arg_genes"].sum().sort_values(ascending=False, kind="stable").index
    keep_cls = list(classes[:top_classes])
    h["class"] = np.where(h["arg_class"].isin(keep_cls), h["arg_class"], "other")
    mat = h[h["host"].isin(hosts_kept)].pivot_table(index="host", columns="class", values="arg_genes", aggfunc="sum",
                                                    fill_value=0)
    mat = mat.reindex(index=hosts_kept, columns=keep_cls + (["other"] if "other" in set(h["class"]) else []),
                      fill_value=0)
    n = len(mat)
    fig, ax = plt.subplots(figsize=fs.figsize("single", n * 3.0 + 28), layout="constrained")
    pal = fs.category_colors(list(mat.columns))
    left = np.zeros(n)
    for c in mat.columns:
        ax.barh(range(n), mat[c].values, left=left, color=pal[c], height=0.75, label=c, lw=0)
        left += mat[c].values
    ax.set_yticks(range(n), mat.index)
    for t in ax.get_yticklabels():
        if t.get_text() != "unclassified" and "(" not in t.get_text():
            t.set_fontstyle("italic")
    ax.set_ylim(n - 0.5, -0.5)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("ARG genes on contigs")
    fig.legend(title="ARG class (DeepARG)", loc="outside lower center", ncol=3, handlelength=1.0, columnspacing=0.8)
    src = mat.rename_axis("host").reset_index().melt(id_vars="host", var_name="arg_class", value_name="arg_genes")
    fw.save(fig, "arg_contig_hosts", "ARG genes on contigs by host",
            f"ARG genes (DeepARG) on assembled contigs per host genus (Metabuli, the contig's lowest rank at genus or "
            f"above in brackets), summed over the assemblies: the {min(top, len(hosts_kept))} hosts with the most "
            f"genes, plus contigs without a classification; colour: the {len(keep_cls)} most common ARG classes, "
            "the rest as 'other'. Physical linkage on one contig, counted in genes, not reads.", src, "contig_args")
