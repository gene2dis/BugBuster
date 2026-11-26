"""argNorm's CARD ARO mapping of ARGs-OAP and DeepARG hits, for bin/results_summary.py (roadmap Phase 27b).

ARGNORM_SARG maps every ARGs-OAP stage two hit (05_arg_prediction/read_level/args_oap/argnorm/
<sample>.args_oap.argnorm.tsv) and ARGNORM_DEEPARG every DeepARG contig hit (contig_level/deeparg/<assembly>/
<assembly>.deeparg.argnorm.tsv) to a CARD ARO, with RGI's cut-off for the mapping (Perfect, Strict, Loose; Manual
for argNorm's own curation) and the drug classes. The summary gives each mapped hit the CARD gene family and
mechanism of the K-MARVEL database's CARD metadata (metadata.tsv = CARD aro_index.tsv, snps.txt), the same gene
key (allele numbers collapsed) and tier as the K-MARVEL / RGI calls. Loose mappings often name a homolog (Phase
27a: DeepARG penA -> mecA): they get no tier and never count as a detection.
"""

from pathlib import Path

import pandas as pd

from summary_args import gene_name, norm_list, tiers_for
from summary_common import die, warn

STRONG = ("Perfect", "Strict", "Manual")
SARG_SUFFIX = ".args_oap.argnorm.tsv"
DEEPARG_SUFFIX = ".deeparg.argnorm.tsv"
NA_LOOSE = "NA_loose_mapping"
NA_NO_CARD = "NA_no_card_metadata"
NA_NOT_IN_CARD = "NA_aro_not_in_card_metadata"
ARGNORM_COLUMNS = ["ARO", "ARO_name", "Cut_Off", "resistance_to_drug_classes_names"]
# ARGs-OAP 4.0.0 blastout.filtered.txt, read by argNorm without its header: columns named 0-17
SARG_ID, SARG_SUBTYPE, SARG_TYPE = "1", "12", "13"
DEEPARG_COLUMNS = ["#ARG", "read_id", "predicted_ARG-class", "best-hit"]


def _read(path, needed):
    """An argNorm output table: its first line is '# argNorm version: <v>', then the header."""
    with open(path) as fh:
        first = fh.readline()
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                     skiprows=1 if first.startswith("# argNorm") else 0)
    if df.empty and not len(df.columns):
        return None
    missing = [c for c in needed if c not in df.columns]
    if missing:
        die(f"{path}: not an argNorm table (lacks {', '.join(missing)})")
    return df


def _classes(text):
    """argNorm's comma-separated drug classes (one per drug it confers resistance to) -> distinct, ';'-joined."""
    seen = []
    for c in str(text).split(","):
        c = c.strip()
        if c and c not in seen:
            seen.append(c)
    return ";".join(seen)


def load_argnorm(root):
    """(SARG hits, samples argNorm ran on; DeepARG hits, assemblies) from one results directory; None, set() for
    a tool without argNorm tables. A table with no hit still counts as run."""
    sdir = root / "05_arg_prediction" / "read_level" / "args_oap" / "argnorm"
    sarg, sran = [], set()
    for p in sorted(sdir.glob(f"*{SARG_SUFFIX}")) if sdir.is_dir() else []:
        sample = p.name[: -len(SARG_SUFFIX)]
        sran.add(sample)
        df = _read(p, ARGNORM_COLUMNS + [SARG_ID, SARG_SUBTYPE, SARG_TYPE])
        if df is None or df.empty:
            continue
        sarg.append(pd.DataFrame({
            "sample": sample, "sarg_type": df[SARG_TYPE], "sarg_subtype": df[SARG_SUBTYPE],
            "sarg_id": df[SARG_ID], "aro_id": df["ARO"], "argnorm_aro_name": df["ARO_name"],
            "cut_off": df["Cut_Off"], "drug_classes": df["resistance_to_drug_classes_names"].map(_classes)}))
    ddir = root / "05_arg_prediction" / "contig_level" / "deeparg"
    deep, dran = [], set()
    paths = (sorted(ddir.glob(f"*/*{DEEPARG_SUFFIX}")) + sorted(ddir.glob(f"coassembly/*/*{DEEPARG_SUFFIX}"))
             if ddir.is_dir() else [])
    for p in paths:
        assembly = p.parent.name
        dran.add(assembly)
        df = _read(p, ARGNORM_COLUMNS + DEEPARG_COLUMNS)
        if df is None or df.empty:
            continue
        deep.append(pd.DataFrame({
            "assembly": assembly, "gene_id": df["read_id"], "contig": df["read_id"].str.rsplit("_", n=1).str[0],
            "arg_gene": df["#ARG"], "arg_class": df["predicted_ARG-class"], "best_hit": df["best-hit"],
            "aro_id": df["ARO"], "argnorm_aro_name": df["ARO_name"], "cut_off": df["Cut_Off"],
            "drug_classes": df["resistance_to_drug_classes_names"].map(_classes)}))
    sarg = pd.concat(sarg, ignore_index=True) if sarg else (pd.DataFrame() if sran else None)
    deep = pd.concat(deep, ignore_index=True) if deep else (pd.DataFrame() if dran else None)
    return sarg, sran, deep, dran


def card_catalog(kmarvel_db):
    """{ARO: (name, AMR gene family, resistance mechanism)} from the K-MARVEL database's CARD metadata
    (metadata.tsv), or None without it."""
    if not kmarvel_db or not (Path(kmarvel_db) / "metadata.tsv").exists():
        return None
    meta = pd.read_csv(Path(kmarvel_db) / "metadata.tsv", sep="\t", dtype=str, keep_default_na=False)
    need = {"ARO Accession", "ARO Name", "AMR Gene Family", "Resistance Mechanism"}
    if not need <= set(meta.columns):
        die(f"{kmarvel_db}/metadata.tsv lacks {', '.join(sorted(need - set(meta.columns)))}")
    meta = meta.drop_duplicates("ARO Accession")
    return {a: (n, norm_list(f), norm_list(m)) for a, n, f, m in
            zip(meta["ARO Accession"], meta["ARO Name"], meta["AMR Gene Family"], meta["Resistance Mechanism"])}


def annotate_card(hits, catalog, families, rules, snps):
    """CARD name, family, gene and tier of each mapped hit (a copy). Loose mappings keep their ARO but get tier
    NA_loose_mapping; without CARD metadata every tier is NA_no_card_metadata and the gene is argNorm's ARO name;
    an ARO the staged CARD release lacks (argNorm maps to CARD 4.0.0) keeps argNorm's name and gets
    NA_aro_not_in_card_metadata."""
    df = hits.copy()
    if df.empty:
        for c in ("aro_name", "gene_family", "gene", "tier", "category", "tier_note"):
            df[c] = pd.Series(dtype=object)
        df["strong"] = pd.Series(dtype=bool)
        return df
    known = df["aro_id"].map(lambda a: catalog.get(a) if catalog else None)
    df["aro_name"] = [k[0] if k else n for k, n in zip(known, df["argnorm_aro_name"])]
    df["gene_family"] = [k[1] if k else "" for k in known]
    mech = [k[2] if k else "" for k in known]
    df["gene"] = [gene_name(n, f, families) if n else "" for n, f in zip(df["aro_name"], df["gene_family"])]
    df["strong"] = df["cut_off"].isin(STRONG) & (df["aro_id"] != "")
    # tier holds 1-4 or an NA code: an object column (pandas 3 would make "" a str column)
    df["tier"] = pd.Series([""] * len(df), index=df.index, dtype=object)
    df["category"], df["tier_note"] = "", ""
    if catalog is None:
        df["tier"] = NA_NO_CARD
        return df
    tierable = df["strong"] & known.notna()
    if tierable.any():
        keys = pd.DataFrame({"aro_id": df.loc[tierable, "aro_id"], "aro_full_name": df.loc[tierable, "aro_name"],
                             "gene_family": df.loc[tierable, "gene_family"], "model": "",
                             "resistance_mechanism": pd.Series(mech, index=df.index)[tierable]})
        t = tiers_for(keys, rules, snps)
        df.loc[tierable, "tier"] = pd.Series([int(x) for x in t["tier"]], index=t.index, dtype=object)
        df.loc[tierable, "category"] = t["category"]
        df.loc[tierable, "tier_note"] = t["tier_note"]
    df.loc[~df["strong"], "tier"] = NA_LOOSE
    df.loc[df["strong"] & known.isna(), "tier"] = NA_NOT_IN_CARD
    return df


def _countable(sarg):
    """The SARG hits that count towards a detection: a strong mapping to an ARO of the CARD metadata (so with a
    CARD gene key and a tier)."""
    return sarg[sarg["strong"] & (sarg["gene"] != "") & sarg["tier"].map(lambda t: isinstance(t, int))]


def sarg_calls(sarg, min_hits):
    """ARGs-OAP's calls for the method comparison, one row per sample x gene (the compare_methods frame):
    detected when at least min_hits hits have a strong (Perfect / Strict / Manual) mapping to the gene; the gene
    takes the lowest tier among them."""
    cols = ["sample", "gene", "tier", "category", "gene_family", "supporting_reads", "gene_representative"]
    if sarg is None or sarg.empty:
        return pd.DataFrame(columns=cols)
    s = _countable(sarg)
    if s.empty:
        return pd.DataFrame(columns=cols)
    s = s.assign(_t=s["tier"].astype(int)).sort_values(["sample", "gene", "_t"], kind="stable")
    agg = s.groupby(["sample", "gene"], sort=False).agg(
        supporting_reads=("aro_id", "size"), gene_family=("gene_family", "first"), tier=("_t", "first"),
        category=("category", "first")).reset_index()
    agg = agg[agg["supporting_reads"] >= min_hits].copy()
    agg["gene_representative"] = True
    return agg[cols].reset_index(drop=True)


def arg_sarg_card_table(sarg, min_hits):
    """arg_sarg_card: one row per sample x SARG subtype x ARO, with the hits, the CARD gene and tier, and whether
    the sample x gene counts as an ARGs-OAP detection in the method comparison."""
    if sarg is None or sarg.empty:
        return None
    keys = ["sample", "sarg_type", "sarg_subtype", "aro_id"]
    out = (sarg.groupby(keys, sort=False)
               .agg(aro_name=("aro_name", "first"), gene=("gene", "first"), gene_family=("gene_family", "first"),
                    cut_off=("cut_off", "first"), read_hits=("sarg_id", "size"), tier=("tier", "first"),
                    category=("category", "first"), drug_classes=("drug_classes", "first"))
               .reset_index())
    n = _countable(sarg).groupby(["sample", "gene"]).size()
    detected = set(n[n >= min_hits].index)
    out["in_comparison"] = [(s, g) in detected and isinstance(t, int) for s, g, t in
                            zip(out["sample"], out["gene"], out["tier"])]
    return out


def arg_deeparg_card_table(parts):
    """arg_deeparg_card: one row per DeepARG hit with its CARD ARO, gene and tier and, when the contig is in a
    bin, the bin. parts: [(batch label, annotated hits, {(assembly, contig): bin}, {(assembly, gene id): bin})]."""
    frames = []
    for label, deep, contig_bin, gene_bin in parts:
        if deep is None or deep.empty:
            continue
        d = deep.copy()
        d.insert(1, "batch", label)
        d["bin"] = [contig_bin.get((a, c)) or gene_bin.get((a, g), "") for a, c, g in
                    zip(d["assembly"], d["contig"], d["gene_id"])]
        frames.append(d)
    if not frames:
        return None
    df = pd.concat(frames, ignore_index=True)
    cols = ["assembly", "batch", "contig", "gene_id", "arg_gene", "arg_class", "best_hit", "aro_id", "aro_name",
            "gene", "gene_family", "cut_off", "tier", "category", "drug_classes", "bin"]
    return df[cols]


def unmapped_warning(sarg, deep):
    """A warning when argNorm left hits without an ARO (argNorm 1.1.0 mapped every SARG 4.0.0 sequence but capO
    and every DeepARG v2 feature, Phase 27a)."""
    for name, df in (("ARGs-OAP", sarg), ("DeepARG", deep)):
        if df is not None and not df.empty:
            n = int((df["aro_id"] == "").sum())
            if n:
                warn(f"argNorm mapped {n} of {len(df)} {name} hits to no ARO; they are left out of the CARD tables")
