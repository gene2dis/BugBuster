"""MAG tables of bin/results_summary.py: one row per refined MAG (CheckM2 quality, MIMAG tier, GTDB-Tk lineage,
Bakta's rRNA / tRNA / CDS counts) and the MAGs' depth per sample.

Published inputs, per assembly <id> (a sample name, coassembly, or a co-assembly group's label; its folder
04_binning/per_sample/<id>, 04_binning/coassembly or 04_binning/coassembly/<group>):
- quality/checkm2/<id>_binette_quality_report.tsv (Binette's CheckM2 report of the refined bins), or with a single
  binner its CheckM2 report <id>_<binner>_quality_report.tsv;
- taxonomy/gtdbtk/<id>_gtdbtk_{bac120,ar53}.tsv;
- 07_functional_annotation/mags/<id>/<id>_<bin>.{tsv,txt} (Bakta, --mag_level_functional; a group's under
  mags/coassembly/<group>/);
- quality/gunc/<id>_gunc_maxCSS_level.tsv (GUNC, --gunc_prediction: one row per refined bin, by its CheckM2 Name);
- coverage/<sample>_bin_depth.tsv (BIN_DEPTH: each MAG's length-weighted mean depth from the jgi depths MetaBAT2
  bins on, one table per sample aligned to the assembly; before roadmap Phase 15 co-assembly only, from
  BEDTOOLS_BIN_DEPTH).

MIMAG (Bowers et al. 2017): high = completeness > 90 %, contamination < 5 %, 23S + 16S + 5S rRNA and >= 18 tRNA
amino acids; medium = completeness >= 50 %, contamination < 10 %; low = the rest. Without Bakta the rRNA / tRNA
part cannot be checked, so a MAG is at most medium (tier_basis checkm2_only). GUNC does not change the tier: its
pass / fail and scores are columns of their own (owner, roadmap Phase 16), added only when GUNC ran.
"""

import re

import numpy as np
import pandas as pd

from summary_common import NA_NOT_RUN, coassembly_subdir, die, load_coassembly_groups, natural_key, warn

RANK_PREFIXES = [("domain", "d__"), ("phylum", "p__"), ("class", "c__"), ("order", "o__"), ("family", "f__"),
                 ("genus", "g__"), ("species", "s__")]
AMINO_ACIDS = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
               "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"}
# GUNC's maxCSS-level table -> mags columns
GUNC_COLUMNS = {"pass.GUNC": "gunc_pass", "clade_separation_score": "clade_separation_score",
                "contamination_portion": "contamination_portion",
                "reference_representation_score": "reference_representation_score"}
CHECKM2_COLUMNS = {"Name", "Completeness", "Contamination", "Completeness_Model_Used", "Coding_Density",
                   "Contig_N50", "Genome_Size", "GC_Content", "Total_Contigs"}


def _assembly_dirs(root, groups=None):
    """{assembly id: its 04_binning folder}; with co-assembly groups (pipeline_info/coassembly_groups.tsv) each
    group's folder is 04_binning/coassembly/<group>."""
    bdir = root / "04_binning"
    out = {}
    groups = load_coassembly_groups(root) if groups is None else groups
    if groups:
        for g in sorted(set(groups.values()), key=natural_key):
            if (bdir / "coassembly" / g).is_dir():
                out[g] = bdir / "coassembly" / g
    elif (bdir / "coassembly").is_dir():
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


def _gunc(gdir, asm):
    """{bin: {gunc column: value}} from the assembly's GUNC table, or None when GUNC did not run. A bin GUNC could
    not score (fewer mapped genes than its --min_mapped_genes) has nan everywhere: blank values."""
    p = gdir / f"{asm}_gunc_maxCSS_level.tsv"
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    if not {"genome", *GUNC_COLUMNS} <= set(df.columns):
        die(f"{p}: GUNC table lacks {sorted({'genome', *GUNC_COLUMNS} - set(df.columns))}")
    out = {}
    for r in df.to_dict("records"):
        passed = {"True": True, "False": False}.get(r["pass.GUNC"])
        scores = {GUNC_COLUMNS[c]: pd.to_numeric(r[c], errors="coerce") for c in GUNC_COLUMNS if c != "pass.GUNC"}
        out[r["genome"]] = {"gunc_pass": passed, **{k: (None if pd.isna(v) else float(v)) for k, v in scores.items()}}
    return out


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
    # header lines start with '#'; pandas' comment= would also cut a field that contains one
    df = pd.read_csv(tsv, sep="\t", header=None, dtype=str, keep_default_na=False,
                     names=["seq", "type", "start", "stop", "strand", "locus", "gene", "product", "dbxrefs"])
    df = df[~df["seq"].str.startswith("#")]
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
    """-> {"mags": DataFrame (one row per refined MAG of every assembly), "depth": DataFrame(assembly, sample, bin,
    mean_depth_x) from the published bin depth tables or None, "bakta": bool} or None when binning did not run."""
    rows, have_bakta, have_gunc, reported = [], False, False, False
    groups = load_coassembly_groups(root)
    for asm, adir in _assembly_dirs(root, groups).items():
        q = _quality_report(asm, adir / "quality" / "checkm2")
        if q is None:
            continue
        reported = True
        binner, qpath = q
        cm = pd.read_csv(qpath, sep="\t", dtype={"Name": str})
        if not CHECKM2_COLUMNS <= set(cm.columns):
            die(f"{qpath}: CheckM2 report lacks {sorted(CHECKM2_COLUMNS - set(cm.columns))}")
        gtdb = _gtdb(adir / "taxonomy" / "gtdbtk", asm)
        gtdb = gtdb.set_index("user_genome") if gtdb is not None else None
        gunc = _gunc(adir / "quality" / "gunc", asm)
        have_gunc |= gunc is not None
        mdir = root / "07_functional_annotation" / "mags" / coassembly_subdir(asm, groups)
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
            if gunc is not None:
                row.update(gunc.get(r.Name, {c: None for c in GUNC_COLUMNS.values()}))
            rows.append(row)
    if not rows:
        # binning ran but no bin passed refinement (header-only reports): an empty table, not a missing stage
        return {"mags": pd.DataFrame(columns=[c for c in MAG_COLUMNS if c != "batch"]), "depth": None,
                "bakta": False, "gunc": have_gunc} if reported else None
    depth = None
    dparts = []
    for asm, adir in _assembly_dirs(root, groups).items():
        cdir = adir / "coverage"
        for p in sorted(cdir.glob("*_bin_depth.tsv")) if cdir.is_dir() else []:
            # BEDTOOLS_BIN_DEPTH before roadmap Phase 11b padded its fields with spaces
            df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
            df = df.apply(lambda c: c.str.strip())
            if not {"sample", "bin_id", "average_cov"} <= set(df.columns):
                die(f"{p}: not a bin depth table (sample, bin_id, Total_contigs, average_cov)")
            dparts.append(pd.DataFrame({"assembly": asm, "sample": df["sample"], "bin": df["bin_id"],
                                        "mean_depth_x": pd.to_numeric(df["average_cov"], errors="coerce")}))
    if dparts:
        depth = pd.concat(dparts, ignore_index=True)
    return {"mags": pd.DataFrame(rows), "depth": depth, "bakta": have_bakta, "gunc": have_gunc, "groups": groups}


MAG_COLUMNS = ["assembly", "batch", "bin", "binner", "completeness_pct", "contamination_pct", "checkm2_model",
               "genome_size_mb", "contigs", "contig_n50_kb", "gc_pct", "coding_density", "quality_tier",
               "tier_basis"] + [r for r, _ in RANK_PREFIXES] + ["gtdb_closest_reference", "rrna_16s", "rrna_23s",
                                                                "rrna_5s", "trna_aa_types", "cds", "pct_hypothetical",
                                                                "crispr_arrays"]


def mag_columns(with_gunc):
    """MAG_COLUMNS, with GUNC's after tier_basis when GUNC ran on any assembly."""
    if not with_gunc:
        return MAG_COLUMNS
    i = MAG_COLUMNS.index("tier_basis") + 1
    return MAG_COLUMNS[:i] + list(GUNC_COLUMNS.values()) + MAG_COLUMNS[i:]


def mags_table(parts):
    """mags: one row per (batch, assembly, refined MAG). parts: [(batch label, load_mags result)]. The GUNC columns
    only when GUNC ran on some assembly; a MAG without a GUNC row (another batch, or an unscored bin) is blank."""
    columns = mag_columns(any(p.get("gunc") for _, p in parts if p is not None))
    frames = [p["mags"].assign(batch=label) for label, p in parts if p is not None]
    if not frames:
        return None
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame(columns=columns)
    df = pd.concat(frames, ignore_index=True).reindex(columns=columns)
    for c in ("contigs", "trna_aa_types", "cds", "crispr_arrays"):
        df[c] = df[c].astype("Int64")
    for c in ("rrna_16s", "rrna_23s", "rrna_5s") + (("gunc_pass",) if "gunc_pass" in columns else ()):
        df[c] = df[c].astype("boolean")
    keys = sorted(range(len(df)), key=lambda i: (df.at[i, "batch"], natural_key(df.at[i, "assembly"]),
                                                 natural_key(df.at[i, "bin"])))
    return df.loc[keys].reset_index(drop=True)


def mag_abundance_frame(parts, contig_tables, sample_bases):
    """Per MAG and sample: mean read depth and the share of the sample's sequenced bases it holds
    (depth x genome size / sample bases x 100).
    - the published bin depth tables (BIN_DEPTH, both modes): the length-weighted mean of the MAG's contig depths
      from the jgi depths MetaBAT2 bins on, per sample aligned to the assembly;
    - a per-sample assembly without them (results of releases before roadmap Phase 15): the length-weighted mean
      of the contig table's coverage over the MAG's contigs, for its own sample (--contig_tax_and_arg);
      NA_not_run without the contig table either.
    parts: [(batch label, load_mags result)]; contig_tables: {batch label: full contig table or None};
    sample_bases: {sample: bases sequenced}. -> DataFrame(assembly, batch, bin, sample, mean_depth_x,
    pct_sample_bases) or None."""
    rows = []
    for label, p in parts:
        if p is None:
            continue
        size = {(r.assembly, r.bin): r.genome_size_mb * 1e6 for r in p["mags"].itertuples(index=False)}
        with_depth = set()
        if p["depth"] is not None:
            for r in p["depth"].itertuples(index=False):
                if (r.assembly, r.bin) in size:
                    rows.append((r.assembly, label, r.bin, r.sample, r.mean_depth_x, size[(r.assembly, r.bin)]))
                    with_depth.add(r.assembly)
        ct = contig_tables.get(label)
        for asm in sorted(set(p["mags"]["assembly"]) - {"coassembly"} - set(p.get("groups", {}).values())
                          - with_depth):
            bins = [b for (a, b) in size if a == asm]
            if ct is None:
                rows += [(asm, label, b, asm, NA_NOT_RUN, size[(asm, b)]) for b in bins]
                continue
            # the contig table has a row per contig and ARG hit: one per contig for the weights
            sub = ct[(ct["Sample"] == asm) & (ct["bin"] != "")].drop_duplicates("contig_name")
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
