"""MAG tables of bin/results_summary.py: one row per refined MAG (CheckM2 quality, MIMAG tier, GTDB-Tk lineage,
Bakta's rRNA / tRNA / CDS counts) and the MAGs' depth per sample.

Published inputs, per assembly <id> (a sample name, or coassembly; its folder 04_binning/per_sample/<id> or
04_binning/coassembly):
- quality/checkm2/<id>_binette_quality_report.tsv (Binette's CheckM2 report of the refined bins), or with a single
  binner its CheckM2 report <id>_<binner>_quality_report.tsv;
- taxonomy/gtdbtk/<id>_gtdbtk_{bac120,ar53}.tsv;
- 07_functional_annotation/mags/<id>/<id>_<bin>.{tsv,txt} (Bakta, --mag_level_functional);
- co-assembly: 04_binning/coassembly/coverage/<sample>_bin_depth.tsv (BEDTOOLS_BIN_DEPTH).

MIMAG (Bowers et al. 2017): high = completeness > 90 %, contamination < 5 %, 23S + 16S + 5S rRNA and >= 18 tRNA
amino acids; medium = completeness >= 50 %, contamination < 10 %; low = the rest. Without Bakta the rRNA / tRNA
part cannot be checked, so a MAG is at most medium (tier_basis checkm2_only).
"""

import re

import numpy as np
import pandas as pd

from summary_common import NA_NOT_RUN, die, natural_key, warn

RANK_PREFIXES = [("domain", "d__"), ("phylum", "p__"), ("class", "c__"), ("order", "o__"), ("family", "f__"),
                 ("genus", "g__"), ("species", "s__")]
AMINO_ACIDS = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
               "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"}
CHECKM2_COLUMNS = {"Name", "Completeness", "Contamination", "Completeness_Model_Used", "Coding_Density",
                   "Contig_N50", "Genome_Size", "GC_Content", "Total_Contigs"}


def _assembly_dirs(root):
    """{assembly id: its 04_binning folder}."""
    bdir = root / "04_binning"
    out = {}
    if (bdir / "coassembly").is_dir():
        out["coassembly"] = bdir / "coassembly"
    for d in sorted((bdir / "per_sample").glob("*")) if (bdir / "per_sample").is_dir() else []:
        if d.is_dir():
            out[d.name] = d
    return out


def _quality_report(asm, qdir):
    """(binner, path) of the assembly's refined-bin CheckM2 report, or None."""
    binette = qdir / f"{asm}_binette_quality_report.tsv"
    if binette.exists():
        return "binette", binette
    reports = sorted(qdir.glob(f"{asm}_*_quality_report.tsv")) if qdir.is_dir() else []
    if len(reports) == 1:
        return reports[0].name[len(asm) + 1:-len("_quality_report.tsv")], reports[0]
    if reports:
        warn(f"{qdir}: several CheckM2 reports and no Binette report; MAGs of {asm} skipped")
    return None


def _gtdb(gdir, asm):
    frames = []
    for dom in ("bac120", "ar53"):
        p = gdir / f"{asm}_gtdbtk_{dom}.tsv"
        if p.exists():
            df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
            if not df.empty:
                if not {"user_genome", "classification"} <= set(df.columns):
                    die(f"{p}: not a GTDB-Tk summary (no user_genome / classification)")
                frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else None


def split_lineage(classification):
    """GTDB 'd__Bacteria;p__...;s__X y' -> {domain..species}; an 'Unclassified ...' call fills domain only."""
    out = {r: "" for r, _ in RANK_PREFIXES}
    parts = {p[:3]: p[3:] for p in str(classification).split(";") if len(p) > 3 and p[1:3] == "__"}
    if not parts and classification:
        out["domain"] = str(classification)
        return out
    for rank, prefix in RANK_PREFIXES:
        out[rank] = parts.get(prefix, "")
    return out


def bakta_features(tsv):
    """Bakta's feature table -> {rrna_16s, rrna_23s, rrna_5s, trna_aa_types}: rRNAs by product name, distinct
    standard amino acids of the tRNAs (product 'tRNA-Met(cat)' -> Met; fMet, SeC and Xxx are not counted)."""
    df = pd.read_csv(tsv, sep="\t", comment="#", header=None, dtype=str, keep_default_na=False,
                     names=["seq", "type", "start", "stop", "strand", "locus", "gene", "product", "dbxrefs"])
    # a truncated rRNA gene ("(3' truncated) 23S ribosomal RNA", at a contig end) counts as present; a
    # pseudo tRNA ("(pseudo) tRNA-Phe(gaa)") does not count
    rrna = {re.sub(r"^\([^)]*\)\s*", "", p) for p in df.loc[df["type"] == "rRNA", "product"]}
    aa = {m.group(1) for p in df.loc[df["type"] == "tRNA", "product"]
          for m in [re.match(r"tRNA-([A-Za-z]{3})", p)] if m and m.group(1) in AMINO_ACIDS}
    return {"rrna_16s": any(p.startswith("16S") for p in rrna), "rrna_23s": any(p.startswith("23S") for p in rrna),
            "rrna_5s": any(p.startswith("5S") for p in rrna), "trna_aa_types": len(aa)}


def bakta_summary(txt):
    """Bakta's .txt summary -> {cds, hypotheticals, crispr_arrays}."""
    vals = {}
    for line in txt.read_text().splitlines():
        key, sep, value = line.partition(":")
        if sep and value.strip().lstrip("-").isdigit():
            vals[key.strip()] = int(value.strip())
    if "CDSs" not in vals:
        die(f"{txt}: no 'CDSs:' line (not a Bakta summary?)")
    return {"cds": vals["CDSs"], "hypotheticals": vals.get("hypotheticals"), "crispr_arrays": vals.get("CRISPR arrays")}


def mimag_tier(completeness, contamination, rrna_all, trna_aa, have_bakta):
    """(quality_tier, tier_basis)."""
    if completeness > 90 and contamination < 5:
        if not have_bakta:
            return "medium", "checkm2_only"
        tier = "high" if rrna_all and trna_aa >= 18 else "medium"
        return tier, "checkm2_rrna_trna"
    basis = "checkm2_rrna_trna" if have_bakta else "checkm2_only"
    if completeness >= 50 and contamination < 10:
        return "medium", basis
    return "low", basis


def load_mags(root):
    """-> {"mags": DataFrame (one row per refined MAG of every assembly), "depth": DataFrame(sample, bin,
    mean_depth_x) of the co-assembly or None, "bakta": bool} or None when binning did not run."""
    rows, have_bakta = [], False
    for asm, adir in _assembly_dirs(root).items():
        q = _quality_report(asm, adir / "quality" / "checkm2")
        if q is None:
            continue
        binner, qpath = q
        cm = pd.read_csv(qpath, sep="\t", dtype={"Name": str})
        if not CHECKM2_COLUMNS <= set(cm.columns):
            die(f"{qpath}: CheckM2 report lacks {sorted(CHECKM2_COLUMNS - set(cm.columns))}")
        gtdb = _gtdb(adir / "taxonomy" / "gtdbtk", asm)
        gtdb = gtdb.set_index("user_genome") if gtdb is not None else None
        mdir = root / "07_functional_annotation" / "mags" / asm
        for r in cm.itertuples(index=False):
            row = {"assembly": asm, "bin": r.Name, "binner": binner,
                   "completeness_pct": float(r.Completeness), "contamination_pct": float(r.Contamination),
                   "checkm2_model": r.Completeness_Model_Used, "genome_size_mb": round(r.Genome_Size / 1e6, 3),
                   "contigs": int(r.Total_Contigs), "contig_n50_kb": round(r.Contig_N50 / 1e3, 3),
                   "gc_pct": round(float(r.GC_Content) * 100, 2), "coding_density": float(r.Coding_Density)}
            g = gtdb.loc[r.Name] if gtdb is not None and r.Name in gtdb.index else None
            if isinstance(g, pd.DataFrame):
                g = g.iloc[0]
            row.update(split_lineage(g["classification"] if g is not None else ""))
            ref = g.get("closest_genome_reference", "") if g is not None else ""
            row["gtdb_closest_reference"] = "" if ref == "N/A" else ref
            tsv, txt = mdir / f"{asm}_{r.Name}.tsv", mdir / f"{asm}_{r.Name}.txt"
            bk = tsv.exists() and txt.exists()
            have_bakta |= bk
            feats = bakta_features(tsv) if bk else {"rrna_16s": None, "rrna_23s": None, "rrna_5s": None,
                                                     "trna_aa_types": None}
            summ = bakta_summary(txt) if bk else {"cds": None, "hypotheticals": None, "crispr_arrays": None}
            row.update(feats)
            row["cds"] = summ["cds"]
            row["pct_hypothetical"] = (round(summ["hypotheticals"] / summ["cds"] * 100, 2)
                                       if bk and summ["hypotheticals"] is not None and summ["cds"] else None)
            row["crispr_arrays"] = summ["crispr_arrays"]
            row["quality_tier"], row["tier_basis"] = mimag_tier(
                row["completeness_pct"], row["contamination_pct"],
                bool(feats["rrna_16s"] and feats["rrna_23s"] and feats["rrna_5s"]), feats["trna_aa_types"] or 0, bk)
            rows.append(row)
    if not rows:
        return None
    depth = None
    cdir = root / "04_binning" / "coassembly" / "coverage"
    dparts = []
    for p in sorted(cdir.glob("*_bin_depth.tsv")) if cdir.is_dir() else []:
        # BEDTOOLS_BIN_DEPTH before roadmap Phase 11b padded its fields with spaces
        df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        df = df.apply(lambda c: c.str.strip())
        if not {"sample", "bin_id", "average_cov"} <= set(df.columns):
            die(f"{p}: not a bin depth table (sample, bin_id, Total_contigs, average_cov)")
        dparts.append(pd.DataFrame({"sample": df["sample"], "bin": df["bin_id"],
                                    "mean_depth_x": pd.to_numeric(df["average_cov"], errors="coerce")}))
    if dparts:
        depth = pd.concat(dparts, ignore_index=True)
    return {"mags": pd.DataFrame(rows), "depth": depth, "bakta": have_bakta}


MAG_COLUMNS = ["assembly", "batch", "bin", "binner", "completeness_pct", "contamination_pct", "checkm2_model",
               "genome_size_mb", "contigs", "contig_n50_kb", "gc_pct", "coding_density", "quality_tier",
               "tier_basis"] + [r for r, _ in RANK_PREFIXES] + ["gtdb_closest_reference", "rrna_16s", "rrna_23s",
                                                                "rrna_5s", "trna_aa_types", "cds", "pct_hypothetical",
                                                                "crispr_arrays"]


def mags_table(parts):
    """mags: one row per (batch, assembly, refined MAG). parts: [(batch label, load_mags result)]."""
    frames = [p["mags"].assign(batch=label) for label, p in parts if p is not None]
    if not frames:
        return None
    df = pd.concat(frames, ignore_index=True)[MAG_COLUMNS]
    for c in ("contigs", "trna_aa_types", "cds", "crispr_arrays"):
        df[c] = df[c].astype("Int64")
    for c in ("rrna_16s", "rrna_23s", "rrna_5s"):
        df[c] = df[c].astype("boolean")
    keys = sorted(range(len(df)), key=lambda i: (df.at[i, "batch"], natural_key(df.at[i, "assembly"]),
                                                 natural_key(df.at[i, "bin"])))
    return df.loc[keys].reset_index(drop=True)


def mag_abundance_frame(parts, contig_tables, sample_bases):
    """Per MAG and sample: mean read depth and the share of the sample's sequenced bases it holds
    (depth x genome size / sample bases x 100).
    - co-assembly: BEDTOOLS_BIN_DEPTH's depth per sample (the mean of the per-contig mean depths, each contig
      weighted equally, as published);
    - per-sample assembly: the length-weighted mean of the contig table's coverage over the MAG's contigs, for
      the MAG's own sample (--contig_tax_and_arg); NA_not_run without the contig table.
    parts: [(batch label, load_mags result)]; contig_tables: {batch label: full contig table or None};
    sample_bases: {sample: bases sequenced}. -> DataFrame(assembly, batch, bin, sample, mean_depth_x,
    pct_sample_bases) or None."""
    rows = []
    for label, p in parts:
        if p is None:
            continue
        size = {(r.assembly, r.bin): r.genome_size_mb * 1e6 for r in p["mags"].itertuples(index=False)}
        if p["depth"] is not None and "coassembly" in set(p["mags"]["assembly"]):
            for r in p["depth"].itertuples(index=False):
                if ("coassembly", r.bin) in size:
                    rows.append(("coassembly", label, r.bin, r.sample, r.mean_depth_x, size[("coassembly", r.bin)]))
        ct = contig_tables.get(label)
        for asm in sorted(set(p["mags"]["assembly"]) - {"coassembly"}):
            bins = [b for (a, b) in size if a == asm]
            if ct is None:
                rows += [(asm, label, b, asm, NA_NOT_RUN, size[(asm, b)]) for b in bins]
                continue
            sub = ct[(ct["Sample"] == asm) & (ct["bin"] != "")]
            length = pd.to_numeric(sub["contig_length"], errors="coerce")
            cov = pd.to_numeric(sub["coverage"], errors="coerce")
            w = (cov * length).groupby(sub["bin"]).sum() / length.groupby(sub["bin"]).sum()
            for b in bins:
                rows.append((asm, label, b, asm, round(float(w[b]), 4) if b in w.index else np.nan, size[(asm, b)]))
    if not rows:
        return None
    df = pd.DataFrame(rows, columns=["assembly", "batch", "bin", "sample", "mean_depth_x", "genome_bp"])

    def share(r):
        if isinstance(r.mean_depth_x, str):
            return r.mean_depth_x
        bases = sample_bases.get(r.sample)
        if bases is None or pd.isna(bases) or not bases or pd.isna(r.mean_depth_x):
            return np.nan
        return round(r.mean_depth_x * r.genome_bp / bases * 100, 4)
    df["pct_sample_bases"] = [share(r) for r in df.itertuples(index=False)]
    return df.drop(columns="genome_bp")
