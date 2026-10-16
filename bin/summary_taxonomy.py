"""Taxonomy loading, merging and watchlist tables for bin/results_summary.py."""

import numpy as np
import pandas as pd

from summary_common import ACCESSION, die, norm_taxon


def load_taxonomy(root):
    tables = root / "02_taxonomy" / "tables"
    otus = sorted(tables.glob("*_otu_table.tsv"))
    if not otus:
        return None
    if len(otus) > 1:
        die(f"{tables} holds more than one *_otu_table.tsv: {[p.name for p in otus]}")
    label = otus[0].name[: -len("_otu_table.tsv")]
    otu = pd.read_csv(otus[0], sep="\t", index_col=0)
    otu.index = otu.index.astype(str)
    tax_path = tables / f"{label}_tax_table.tsv"
    if not tax_path.exists():
        die(f"{tax_path} is missing")
    tax = pd.read_csv(tax_path, sep="\t", index_col=0, dtype=str, keep_default_na=False)
    tax.index = tax.index.astype(str)
    counts_path = tables / f"{label}_read_counts.tsv"
    counts = pd.read_csv(counts_path, sep="\t", index_col=0) if counts_path.exists() else None
    if counts is not None:
        counts.index = counts.index.astype(str)
    sylph = load_sylph_details(root, label)
    return {"label": label, "otu": otu, "tax": tax, "counts": counts, "sylph": sylph}


def load_sylph_details(root, label):
    """Per-genome ANI and coverage from the raw sylph profile TSVs (sylph runs only)."""
    sylph_dir = root / "02_taxonomy" / "sylph"
    if not label.startswith("sylph_") or not sylph_dir.is_dir():
        return None
    suffix = f"_{label}.tsv"
    frames = []
    for path in sorted(sylph_dir.glob(f"*{suffix}")):
        df = pd.read_csv(path, sep="\t", dtype=str)
        if df.empty:
            continue
        acc = df["Genome_file"].str.extract(ACCESSION, expand=False)
        fallback = df["Genome_file"].str.replace(r".*/", "", regex=True).str.replace(
            r"(_genomic)?\.(fna|fa|fasta)(\.gz)?$", "", regex=True)
        ci = df.get("ANI_5-95_percentile", pd.Series("NA-NA", index=df.index)).fillna("NA-NA")
        ci_low = pd.to_numeric(ci.str.split("-").str[0], errors="coerce")
        frames.append(pd.DataFrame({
            "sample": path.name[: -len(suffix)],
            "taxon_id": acc.fillna(fallback),
            "sequence_abundance": pd.to_numeric(df["Sequence_abundance"], errors="coerce"),
            "adjusted_ani": pd.to_numeric(df["Adjusted_ANI"], errors="coerce"),
            "ani_5_95": ci.where(ci != "NA-NA", ""),
            "true_cov": pd.to_numeric(df["True_cov"], errors="coerce"),
            "eff_lambda": df["Eff_lambda"],
            "species_call_uncertain": ci_low < 95,
        }))
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def load_watchlist(path):
    wl = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    if not {"taxon", "rank"} <= set(wl.columns):
        die(f"{path} needs columns taxon, rank")
    if "note" not in wl.columns:
        wl["note"] = ""
    return wl


def build_taxonomy(batches, samples, order):
    tax_batches = [b["taxonomy"] for b in batches if b["taxonomy"] is not None]
    if not tax_batches:
        return None
    labels = {t["label"] for t in tax_batches}
    if len(labels) > 1:
        die(f"results directories use different taxonomy profilers/databases: {sorted(labels)}")
    otu = pd.concat([t["otu"] for t in tax_batches], axis=1).fillna(0)
    tax = pd.concat([t["tax"] for t in tax_batches])
    tax = tax[~tax.index.duplicated()].reindex(otu.index).fillna("")
    counts = None
    if all(t["counts"] is not None for t in tax_batches):
        counts = pd.concat([t["counts"] for t in tax_batches], axis=1).reindex(otu.index).fillna(0)
        counts = counts.round().astype("int64")  # the outer join made them float
    sylph = [t["sylph"] for t in tax_batches if t["sylph"] is not None]
    sylph = pd.concat(sylph, ignore_index=True) if sylph else None
    totals = otu.sum(axis=0).replace(0, np.nan)
    rel = otu.div(totals, axis=1) * 100
    cols = [s for s in order if s in rel.columns]
    rel = rel.reindex(columns=cols)
    ranks = list(tax.columns)

    long = rel.stack().rename("relative_abundance").reset_index()
    long.columns = ["taxon_id", "sample", "relative_abundance"]
    long = long[long["relative_abundance"] > 0]
    if counts is not None:
        cstack = counts.stack().rename("est_read_pairs").reset_index()
        cstack.columns = ["taxon_id", "sample", "est_read_pairs"]
        long = long.merge(cstack, on=["taxon_id", "sample"], how="left")
    if sylph is not None:
        long = long.merge(sylph, on=["sample", "taxon_id"], how="left")
    long = long.merge(tax, left_on="taxon_id", right_index=True, how="left")

    wide = {}
    for rank in ranks:
        names = tax[rank].replace("", np.nan).fillna(f"unassigned {rank.lower()}")
        agg = rel.groupby(names).sum()
        agg = agg.loc[agg.mean(axis=1).sort_values(ascending=False).index]
        agg.index.name = rank
        wide[rank] = agg
    return {"label": labels.pop(), "rel": rel, "tax": tax, "long": long, "wide": wide, "ranks": ranks,
            "sylph": sylph is not None}


def species_rank(tx):
    for cand in ("Species", "species"):
        if cand in tx["wide"]:
            return cand
    return tx["ranks"][-1]


def genus_rank(tx):
    for cand in ("Genus", "genus"):
        if cand in tx["wide"]:
            return cand
    return tx["ranks"][-2] if len(tx["ranks"]) > 1 else tx["ranks"][-1]


def prevalence_table(wide, samples, group_by):
    present = wide > 0
    out = pd.DataFrame({
        "n_samples_detected": present.sum(axis=1),
        "pct_samples_detected": present.mean(axis=1) * 100,
        "mean_abundance": wide.mean(axis=1),
        "median_abundance_when_detected": wide.where(present).median(axis=1),
        "max_abundance": wide.max(axis=1),
    })
    if group_by:
        groups = samples.loc[wide.columns, group_by]
        out.insert(2, f"n_{group_by}_detected", present.T.groupby(groups.values).any().sum(axis=0))
    out = out[out["n_samples_detected"] > 0].sort_values(["n_samples_detected", "mean_abundance"],
                                                          ascending=False)
    return out


def watchlist_table(tx, watch):
    sp, ge = species_rank(tx), genus_rank(tx)
    long = tx["long"]
    rows, long_rows = [], []
    for w in watch.itertuples(index=False):
        rank_col = sp if w.rank.lower() == "species" else ge
        hit = long[long[rank_col].map(norm_taxon) == norm_taxon(w.taxon)]
        per_sample = hit.groupby("sample")["relative_abundance"].sum()
        names = sorted(set(hit[rank_col]))
        row = {"taxon": w.taxon, "rank": w.rank, "note": w.note, "matched_names": "; ".join(names),
               "n_samples_detected": int((per_sample > 0).sum())}
        for s in tx["rel"].columns:
            row[s] = float(per_sample.get(s, 0.0))
        rows.append(row)
        for s, ab in per_sample.items():
            sub = hit[hit["sample"] == s]
            uncertain = bool(sub["species_call_uncertain"].fillna(False).all()) if "species_call_uncertain" in sub else False
            long_rows.append({"taxon": w.taxon, "sample": s, "relative_abundance": ab, "uncertain": uncertain})
    return pd.DataFrame(rows), pd.DataFrame(long_rows)


def bray_curtis_pcoa(rel, n_axes=3):
    """Relative abundance (taxa x samples) -> (coordinates: samples x PC1..PCn, % variance of each axis).

    Bray-Curtis dissimilarity between samples (on proportions), then classical PCoA (Gower double
    centring). Bray-Curtis is not Euclidean, so some eigenvalues are negative; the % variance of an
    axis is its eigenvalue over the sum of the positive ones.
    """
    x = rel.fillna(0).T.astype(float)
    x = x[x.sum(axis=1) > 0]
    if len(x) < 3:
        return None, None
    p = x.values / x.values.sum(axis=1, keepdims=True)
    n = len(p)
    d = np.zeros((n, n))
    for i in range(n):
        d[i] = np.abs(p - p[i]).sum(axis=1) / (p + p[i]).sum(axis=1)
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d ** 2) @ j
    vals, vecs = np.linalg.eigh(b)
    idx = np.argsort(vals)[::-1]
    vals, vecs = vals[idx], vecs[:, idx]
    pos = vals[vals > 1e-12]
    k = min(n_axes, len(pos))
    coords = pd.DataFrame(vecs[:, :k] * np.sqrt(vals[:k]), index=x.index,
                          columns=[f"PC{i + 1}" for i in range(k)])
    pct = list(vals[:k] / pos.sum() * 100)
    return coords, pct
