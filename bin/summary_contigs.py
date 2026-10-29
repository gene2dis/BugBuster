"""Contig-level ARGs for bin/results_summary.py: ARG genes on geNomad's plasmid and virus contigs and inside its
proviruses, from ARG_CONTIG_LEVEL_REPORT's Contig_tax_and_arg_prediction.tsv (--contig_tax_and_arg with
--genomad_prediction)."""

import pandas as pd

CONTIG_TABLE = ("05_arg_prediction", "contig_level", "summary", "Contig_tax_and_arg_prediction.tsv")
NEEDED = {"Sample", "contig_name", "ARG_gene", "predicted_ARG_class"}
GENOMAD = {"mobile_class", "provirus"}


def load_contig_args(root):
    """The contig table's ARG rows (one per ARG gene), or None when the stage did not run. The geNomad columns
    (mobile_class, provirus) are there only when --genomad_prediction ran too."""
    path = root.joinpath(*CONTIG_TABLE)
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = NEEDED - set(df.columns)
    if missing:
        raise SystemExit(f"ERROR: {path} lacks columns {sorted(missing)}")
    return df[df["ARG_gene"] != ""].reset_index(drop=True)


def has_genomad(df):
    """The contig table carries geNomad's columns."""
    return df is not None and GENOMAD <= set(df.columns)


def contig_mobile_table(parts):
    """arg_contig_mobile: per assembly (sample, or coassembly) and DeepARG class, the ARG genes and how many sit
    on plasmid contigs, on virus contigs and inside proviruses. None when no batch has the geNomad columns."""
    parts = [p for p in parts if has_genomad(p) and not p.empty]
    if not parts:
        return None
    df = pd.concat(parts, ignore_index=True)
    df["on_plasmid"] = df["mobile_class"] == "plasmid"
    df["on_virus"] = df["mobile_class"] == "virus"
    df["in_provirus"] = df["provirus"] != ""
    df["mobile"] = df["on_plasmid"] | df["on_virus"] | df["in_provirus"]
    out = (df.groupby(["Sample", "predicted_ARG_class"])
             .agg(arg_genes=("ARG_gene", "size"), arg_genes_on_plasmid=("on_plasmid", "sum"),
                  arg_genes_on_virus=("on_virus", "sum"), arg_genes_in_provirus=("in_provirus", "sum"),
                  mobile=("mobile", "sum"))
             .reset_index()
             .rename(columns={"Sample": "assembly", "predicted_ARG_class": "arg_class"}))
    out["pct_mobile"] = (out.pop("mobile") / out["arg_genes"] * 100).round(2)
    return out.sort_values(["assembly", "arg_genes"], ascending=[True, False], kind="stable").reset_index(drop=True)
