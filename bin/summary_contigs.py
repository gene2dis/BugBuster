"""Contig-level ARGs for bin/results_summary.py, from ARG_CONTIG_LEVEL_REPORT's
Contig_tax_and_arg_prediction.tsv (--contig_tax_and_arg; geNomad columns with --genomad_prediction): ARG genes on
geNomad's plasmid and virus contigs and inside its proviruses, every ARG hit with the taxon and bin of its contig,
ARG genes per host taxon, and (with --contig_level_functional) each ARG gene's abundance per sample from
DEEPARG_CONTIGS' hits and AGGREGATE_FUNCTIONS' gene_abundance.tsv.

This is linkage: the ARG and the taxon are on the same assembled contig. Counts are genes on contigs, not read
abundances."""

import pandas as pd

from summary_common import die

CONTIG_TABLE = ("05_arg_prediction", "contig_level", "summary", "Contig_tax_and_arg_prediction.tsv")
NEEDED = {"Sample", "contig_name", "ARG_gene", "predicted_ARG_class"}
GENOMAD = {"mobile_class", "provirus"}
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
MAPPING_COLUMNS = ["#ARG", "read_id", "predicted_ARG-class"]
GENE_ABUNDANCE_COLUMNS = ["sample_id", "gene_id", "assembly_mode", "count", "length_bp", "tpm", "cpge"]


def load_contig_table(root):
    """Every row of the contig table (one per contig and ARG hit; contigs without an ARG have one row with an empty
    ARG_gene), or None when the stage did not run."""
    path = root.joinpath(*CONTIG_TABLE)
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = NEEDED - set(df.columns)
    if missing:
        raise SystemExit(f"ERROR: {path} lacks columns {sorted(missing)}")
    return df


def arg_rows(table):
    """The contig table's ARG rows (one per ARG gene), or None."""
    if table is None:
        return None
    return table[table["ARG_gene"] != ""].reset_index(drop=True)


def load_contig_args(root):
    """The contig table's ARG rows (one per ARG gene), or None when the stage did not run. The geNomad columns
    (mobile_class, provirus) are there only when --genomad_prediction ran too."""
    return arg_rows(load_contig_table(root))


def has_genomad(df):
    """The contig table carries geNomad's columns."""
    return df is not None and GENOMAD <= set(df.columns)


def contig_mobile_table(parts):
    """arg_contig_mobile: per batch, assembly (sample, or coassembly) and DeepARG class, the ARG genes and how many
    sit on plasmid contigs, on virus contigs and inside proviruses. parts: [(batch label, contig ARG rows)]. The
    batch is part of the key: every co-assembly is called 'coassembly'. None when no batch has the geNomad
    columns."""
    parts = [p.assign(batch=label) for label, p in parts if has_genomad(p) and not p.empty]
    if not parts:
        return None
    df = pd.concat(parts, ignore_index=True)
    df["on_plasmid"] = df["mobile_class"] == "plasmid"
    df["on_virus"] = df["mobile_class"] == "virus"
    df["in_provirus"] = df["provirus"] != ""
    df["mobile"] = df["on_plasmid"] | df["on_virus"] | df["in_provirus"]
    out = (df.groupby(["batch", "Sample", "predicted_ARG_class"])
             .agg(arg_genes=("ARG_gene", "size"), arg_genes_on_plasmid=("on_plasmid", "sum"),
                  arg_genes_on_virus=("on_virus", "sum"), arg_genes_in_provirus=("in_provirus", "sum"),
                  mobile=("mobile", "sum"))
             .reset_index()
             .rename(columns={"Sample": "assembly", "predicted_ARG_class": "arg_class"}))
    out["pct_mobile"] = (out.pop("mobile") / out["arg_genes"] * 100).round(2)
    out = out[["assembly", "batch"] + [c for c in out.columns if c not in ("assembly", "batch")]]
    return out.sort_values(["batch", "assembly", "arg_genes"], ascending=[True, True, False],
                           kind="stable").reset_index(drop=True)


def host(row):
    """(rank, taxon) of the lowest rank Metabuli filled for the contig, or ('unclassified', '')."""
    for rank in reversed(RANKS):
        if row.get(rank, ""):
            return rank, row[rank]
    return "unclassified", ""


def _species(gtdb):
    return next((p[3:] for p in str(gtdb).split(";") if p.startswith("s__")), "")


def arg_contig_long_table(parts, mags=None):
    """arg_contig_long: one row per DeepARG hit on a contig of each (batch, assembly), with the contig's length,
    depth, Metabuli taxon (lowest rank reached as host_rank / host_taxon), bin, the bin's MIMAG tier (from mags)
    and GTDB species, and geNomad's class. parts: [(batch label, contig ARG rows)]."""
    frames = []
    for label, df in parts:
        if df is None or df.empty:
            continue
        hosts = [host(r) for r in df.to_dict("records")]
        out = pd.DataFrame({
            "assembly": df["Sample"], "batch": label, "contig": df["contig_name"],
            "contig_length_bp": pd.to_numeric(df.get("contig_length"), errors="coerce"),
            "contig_depth_x": pd.to_numeric(df.get("coverage"), errors="coerce"),
            "arg_gene": df["ARG_gene"], "arg_class": df["predicted_ARG_class"],
            "deeparg_probability": pd.to_numeric(df.get("probability"), errors="coerce"),
            "identity_pct": pd.to_numeric(df.get("identity"), errors="coerce"),
            "host_rank": [h[0] for h in hosts], "host_taxon": [h[1] for h in hosts],
            "phylum": df.get("phylum", ""), "genus": df.get("genus", ""), "species": df.get("species", ""),
            "taxonomy_score": pd.to_numeric(df.get("taxonomy_score"), errors="coerce"),
            "bin": df.get("bin", ""),
            "bin_gtdb_species": df["bin_gtdb_classification"].map(_species) if "bin_gtdb_classification" in df else "",
            "mobile_class": df.get("mobile_class", ""),
        })
        frames.append(out)
    if not frames:
        return None
    df = pd.concat(frames, ignore_index=True)
    tier = {}
    if mags is not None:
        tier = {(r.batch, r.assembly, r.bin): r.quality_tier for r in mags.itertuples(index=False)}
    df.insert(df.columns.get_loc("bin") + 1, "mag_tier",
              [tier.get((b, a, x), "") if x else "" for b, a, x in zip(df["batch"], df["assembly"], df["bin"])])
    df["contig_length_bp"] = df["contig_length_bp"].astype("Int64")
    df = df.sort_values(["batch", "assembly", "arg_class", "arg_gene", "contig"], kind="stable")
    return df.reset_index(drop=True)


def arg_contig_hosts_table(long):
    """arg_contig_hosts: per (batch, assembly, ARG class, host): ARG genes, the contigs carrying them, those on a
    binned contig and those on a geNomad plasmid or virus contig. The host is the contig's genus, or the lowest rank
    above it Metabuli reached (species calls are counted at their genus)."""
    if long is None or long.empty:
        return None
    df = long.copy()
    sp = (df["host_rank"] == "species") & (df["genus"] != "")
    df.loc[sp, "host_rank"] = "genus"
    df.loc[sp, "host_taxon"] = df.loc[sp, "genus"]
    df["on_mag"] = df["bin"] != ""
    df["on_mobile"] = df["mobile_class"].isin(["plasmid", "virus"])
    out = (df.groupby(["batch", "assembly", "arg_class", "host_rank", "host_taxon"])
             .agg(arg_genes=("arg_gene", "size"), arg_contigs=("contig", "nunique"), on_mag=("on_mag", "sum"),
                  on_mobile=("on_mobile", "sum"))
             .reset_index())
    out = out[["assembly", "batch", "arg_class", "host_rank", "host_taxon", "arg_genes", "arg_contigs", "on_mag",
               "on_mobile"]]
    return out.sort_values(["batch", "assembly", "arg_class", "arg_genes"], ascending=[True, True, True, False],
                           kind="stable").reset_index(drop=True)


def load_arg_gene_abundance(root):
    """DEEPARG_CONTIGS' hits per assembly (05_arg_prediction/contig_level/deeparg/<id>/*.mapping.ARG, read_id =
    the Pyrodigal gene id) joined to gene_abundance.tsv (07_functional_annotation/summary/), keyed on the assembly:
    gene ids repeat across per-sample assemblies (MEGAHIT names every assembly's contigs k141_N), so a per-sample
    assembly's genes are matched to its own sample only, a co-assembly's to every sample. None without both."""
    ddir = root / "05_arg_prediction" / "contig_level" / "deeparg"
    gpath = root / "07_functional_annotation" / "summary" / "gene_abundance.tsv"
    maps = sorted(ddir.glob("*/*.mapping.ARG")) if ddir.is_dir() else []
    if not maps or not gpath.exists():
        return None
    hits = []
    for p in maps:
        df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        if not set(MAPPING_COLUMNS) <= set(df.columns):
            die(f"{p}: not a DeepARG mapping.ARG table (lacks {sorted(set(MAPPING_COLUMNS) - set(df.columns))})")
        hits.append(pd.DataFrame({"assembly": p.parent.name, "gene_id": df["read_id"], "arg_gene": df["#ARG"],
                                  "arg_class": df["predicted_ARG-class"]}))
    hits = pd.concat(hits, ignore_index=True)
    if hits.empty:
        return None
    wanted = set(hits["gene_id"])
    head = pd.read_csv(gpath, sep="\t", nrows=0).columns.tolist()
    if head != GENE_ABUNDANCE_COLUMNS:
        die(f"{gpath}: unexpected header {head} (expected {GENE_ABUNDANCE_COLUMNS})")
    chunks = [c[c["gene_id"].isin(wanted)]
              for c in pd.read_csv(gpath, sep="\t", dtype={"sample_id": str, "gene_id": str}, chunksize=500_000)]
    ab = pd.concat(chunks, ignore_index=True)
    coassembly = (ab["assembly_mode"] == "coassembly").any() if not ab.empty else "coassembly" in set(hits["assembly"])
    ab["assembly"] = "coassembly" if coassembly else ab["sample_id"]
    out = hits.merge(ab, on=["assembly", "gene_id"], how="inner")
    out = out.rename(columns={"sample_id": "sample"})
    out["contig"] = out["gene_id"].str.rsplit("_", n=1).str[0]
    return out[["sample", "assembly", "gene_id", "contig", "arg_gene", "arg_class", "count", "tpm", "cpge"]]
