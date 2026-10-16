#!/usr/bin/env python3
"""
BugBuster results summary: consolidated, cross-batch tables, figures and an HTML report.

Usage:
  results_summary.py --results RUN1/results [RUN2/results ...] --outdir summary
                     [--metadata samples.tsv --group-by subject --order-by visit]
                     [--kmarvel-db DIR] [--card-prevalence card_prevalence.txt]

Reads only the files a BugBuster run publishes, so it runs on finished results without
rerunning any tool, and merges any number of results directories (one per batch):
  01_quality_control/summary/Reads_report.csv     read counts per step (single reads)
  01_quality_control/fastp/*/*.fastp.json         mean read length
  02_taxonomy/tables/<label>_{otu,tax}_table.tsv   any taxonomic profiler (+ _read_counts.tsv)
  02_taxonomy/sylph/*_sylph_<db>.tsv              sylph ANI / coverage (optional)
  05_arg_prediction/read_level/kmarvel/*/*_kmarvel.tsv   K-MARVEL ARG calls
  pipeline_info/software_versions.yml             tool and database versions
Each section is skipped when its inputs are absent.

Units: Reads_report.csv counts single reads (both mates); every count here is in read
pairs (= reads / 2), the unit K-MARVEL and sylph use. ARG abundance is normalised to the
sample's clean read pairs: RPM = supporting read pairs / clean pairs x 1e6, and
RPKM = RPM / (reference gene length in kb), the length being the CARD reference
protein x 3. Samples below --min-clean-pairs are flagged (never dropped): an ARG or
taxon not found there is reported as NA_low_depth, not as absent.

ARGs are assigned a reporting tier from a user-editable rule table (assets/arg_priority_tiers.tsv):
1 clinical priority, 2 other acquired, 3 intrinsic/efflux/regulator, 4 mutation dependent.
Tier 4 calls are checked against CARD snps.txt (from --kmarvel-db): K-MARVEL's
mutation_data lists every difference from the reference, so a variant-model hit is a
resistance call only when a known resistance mutation is among them.

Outputs (in --outdir): tables/*.tsv, results_summary.xlsx (if openpyxl or xlsxwriter
is installed), figures/*.{png,pdf,svg} and results_summary.html (self-contained, works offline).
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# the summary_* modules sit next to this script; put its directory on the path
# also when the script is loaded by file path (importlib) rather than run, and
# write no __pycache__ into the pipeline's bin/ (read-only in a container image)
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from summary_common import NA_NOT_RUN, die, natural_key, norm_taxon, strip_gtdb, warn  # noqa: E402
from summary_taxonomy import (build_taxonomy, load_taxonomy, load_watchlist, prevalence_table,  # noqa: E402
                              species_rank, watchlist_table)
from summary_args import (annotate_args, arg_status, arg_tables, load_kmarvel, load_prevalence,  # noqa: E402
                          load_rules, load_snps)
from summary_figures import FigureWriter, fig_args, fig_qc, fig_taxonomy, fig_watchlist  # noqa: E402
from summary_html import write_excel, write_html, write_tables  # noqa: E402


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_ASSETS = SCRIPT_DIR.parent / "assets"


def batch_label(path):
    p = Path(path).resolve()
    return p.parent.name if p.name == "results" else p.name


def load_reads_report(root):
    path = root / "01_quality_control" / "summary" / "Reads_report.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path, dtype={"Id": str})
    if "Id" not in df.columns:
        die(f"{path} has no Id column")
    steps = [
        c for c in df.columns
        if c != "Id" and "singleton" not in c.lower()
        and c not in ("Classified", "Unclassified")
        and pd.api.types.is_numeric_dtype(df[c])
    ]
    if not steps:
        die(f"{path} has no read-count columns")
    raw = "Raw reads" if "Raw reads" in steps else steps[0]
    qc = "Fastp" if "Fastp" in steps else raw
    clean = steps[-1]
    out = pd.DataFrame({
        "sample": df["Id"],
        "raw_pairs": df[raw] // 2,
        "qc_pass_pairs": df[qc] // 2,
        "clean_pairs": df[clean] // 2,
    })
    if "Classified" in df.columns:
        classified = pd.to_numeric(df["Classified"], errors="coerce")
        # TAXONOMY_REPORT writes kraken2's as a percentage (0-100) and sylph's /
        # sourmash's as a fraction; the profiler's DB column tells them apart
        # (a value test would misread a batch where every sample is <= 1% classified)
        if "Kraken DB" in df.columns:
            classified = classified / 100
        out["classified_fraction"] = classified
    return out


def load_fastp_lengths(root):
    rows = []
    for path in sorted((root / "01_quality_control" / "fastp").glob("*/*.fastp.json")):
        try:
            after = json.loads(path.read_text())["summary"]["after_filtering"]
        except (KeyError, ValueError) as err:
            warn(f"cannot read {path}: {err}")
            continue
        lengths = [after.get("read1_mean_length"), after.get("read2_mean_length")]
        lengths = [x for x in lengths if x]
        rows.append({
            "sample": path.name[: -len(".fastp.json")],
            "mean_read_length": float(np.mean(lengths)) if lengths else np.nan,
        })
    return pd.DataFrame(rows, columns=["sample", "mean_read_length"])


def load_versions(root):
    path = root / "pipeline_info" / "software_versions.yml"
    return path.read_text() if path.exists() else ""


def load_metadata(path):
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    meta = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    key = "sample" if "sample" in meta.columns else meta.columns[0]
    meta = meta.rename(columns={key: "sample"})
    if meta["sample"].duplicated().any():
        die(f"{path}: duplicate sample names: {sorted(meta.loc[meta['sample'].duplicated(), 'sample'])[:5]}")
    reserved = {"batch"}
    clash = reserved & set(meta.columns)
    if clash:
        meta = meta.rename(columns={c: f"meta_{c}" for c in clash})
    return meta


def sample_order(samples, group_by, order_by):
    keys = []
    for s in samples.index:
        g = samples.at[s, group_by] if group_by else ""
        o = samples.at[s, order_by] if order_by else ""
        keys.append((natural_key(g), natural_key(o), natural_key(s), s))
    return [k[-1] for k in sorted(keys)]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", nargs="+", required=True, help="BugBuster results directories (one per batch)")
    p.add_argument("--batch-labels", nargs="+", metavar="LABEL",
                   help="batch labels, one per --results directory in the same order "
                        "(default: each directory's name, or its parent's when it is called 'results')")
    p.add_argument("--outdir", required=True, help="output directory")
    p.add_argument("--metadata", help="TSV/CSV with a 'sample' column (or sample names in the first column)")
    p.add_argument("--group-by", help="metadata column that groups samples (e.g. subject, site)")
    p.add_argument("--order-by", help="metadata column that orders samples within a group (e.g. visit, time)")
    p.add_argument("--min-clean-pairs", type=int, default=1_000_000,
                   help="flag samples below this many clean read pairs as low depth (default 1,000,000)")
    p.add_argument("--kmarvel-db", help="K-MARVEL database directory (for CARD snps.txt)")
    p.add_argument("--card-prevalence", help="WildCARD card_prevalence.txt, for candidate ARG hosts")
    p.add_argument("--arg-priority", default=str(DEFAULT_ASSETS / "arg_priority_tiers.tsv"),
                   help="ARG tier rule table (default: assets/arg_priority_tiers.tsv)")
    p.add_argument("--watchlist", default=str(DEFAULT_ASSETS / "summary_watchlist.tsv"),
                   help="watchlist taxa (default: assets/summary_watchlist.tsv)")
    p.add_argument("--top-n", type=int, default=40, help="maximum rows in the taxonomy and ARG heatmaps (default 40)")
    p.add_argument("--title", default="BugBuster results summary", help="HTML report title")
    p.add_argument("--no-excel", action="store_true", help="do not write the Excel workbook")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.order_by and not args.group_by:
        warn("--order-by without --group-by orders all samples as one series")

    if args.batch_labels is not None:
        if len(args.batch_labels) != len(args.results):
            die(f"--batch-labels has {len(args.batch_labels)} labels for {len(args.results)} --results directories")
        if any(not lab.strip() for lab in args.batch_labels):
            die("--batch-labels: a label is empty")
        dups = sorted({lab for lab in args.batch_labels if args.batch_labels.count(lab) > 1})
        if dups:
            die(f"--batch-labels: {', '.join(dups)} given more than once; batch labels must be unique")

    batches, seen, labels = [], {}, set()
    for i, r in enumerate(args.results):
        root = Path(r)
        if not root.is_dir():
            die(f"{root} is not a directory")
        if args.batch_labels is not None:
            label = args.batch_labels[i]
        else:
            label = batch_label(root)
            if label in labels:
                label = f"{label}_{i + 1}"
        labels.add(label)
        reads = load_reads_report(root)
        fastp = load_fastp_lengths(root)
        taxonomy = load_taxonomy(root)
        kmarvel, kran = load_kmarvel(root)
        names = set()
        if reads is not None:
            names |= set(reads["sample"])
        if taxonomy is not None:
            names |= set(taxonomy["otu"].columns)
        names |= kran
        if not names:
            die(f"{root}: no BugBuster outputs found (expected 01_quality_control/, 02_taxonomy/ or 05_arg_prediction/)")
        for s in names:
            if s in seen:
                die(f"sample {s} is in both {seen[s]} and {root}; sample names must be unique across results directories")
            seen[s] = str(root)
        batches.append({"label": label, "root": root, "reads": reads, "fastp": fastp, "taxonomy": taxonomy,
                        "kmarvel": kmarvel, "kmarvel_ran": kran, "samples": names, "versions": load_versions(root)})
        print(f"{root}: {len(names)} samples (batch '{label}')")

    samples = pd.DataFrame(
        [{"sample": s, "batch": b["label"]} for b in batches for s in sorted(b["samples"])]).set_index("sample")
    meta_cols, key_meta = [], []
    if args.metadata:
        meta = load_metadata(args.metadata).set_index("sample")
        unknown = sorted(set(meta.index) - set(samples.index))
        missing = sorted(set(samples.index) - set(meta.index))
        if unknown:
            warn(f"{len(unknown)} metadata samples are not in the results (e.g. {unknown[:3]})")
        if missing:
            warn(f"{len(missing)} samples have no metadata (e.g. {missing[:3]})")
        meta_cols = list(meta.columns)
        samples = samples.join(meta, how="left")
        samples[meta_cols] = samples[meta_cols].fillna("")
    for col in (args.group_by, args.order_by):
        if col and col not in samples.columns:
            die(f"--group-by/--order-by column '{col}' is not in the metadata ({meta_cols or 'no --metadata given'})")
    key_meta = [c for c in (args.group_by, args.order_by) if c]
    order = sample_order(samples, args.group_by, args.order_by)
    samples = samples.loc[order]
    groups = [samples.at[s, args.group_by] if args.group_by else "" for s in order]

    # ---- QC
    reads = [b["reads"] for b in batches if b["reads"] is not None]
    if reads:
        r = pd.concat(reads).set_index("sample")
        samples = samples.join(r, how="left")
        samples["pct_removed_qc"] = (1 - samples["qc_pass_pairs"] / samples["raw_pairs"]) * 100
        samples["pct_removed_decontam"] = (samples["qc_pass_pairs"] - samples["clean_pairs"]) / samples["raw_pairs"] * 100
        samples["pct_clean"] = samples["clean_pairs"] / samples["raw_pairs"] * 100
        samples["low_depth"] = samples["clean_pairs"] < args.min_clean_pairs
    else:
        warn("no Reads_report.csv: depth flags and ARG normalisation are unavailable")
    lengths = pd.concat([b["fastp"] for b in batches])
    if not lengths.empty:
        samples = samples.join(lengths.set_index("sample"), how="left")
    low = set(samples.index[samples["low_depth"].fillna(False).astype(bool)]) if "low_depth" in samples else set()

    tables = {}
    figs = FigureWriter(outdir)
    summary = {"n_samples": len(order), "n_batches": len(batches), "n_low_depth": len(low)}

    # ---- taxonomy
    tx = build_taxonomy(batches, samples, order)
    species_by_sample = {}
    if tx is not None:
        sp = species_rank(tx)
        spw = tx["wide"][sp].reindex(columns=order).fillna(0)
        p = spw / 100
        with np.errstate(divide="ignore", invalid="ignore"):
            shannon = -(p * np.log(p.where(p > 0))).sum(axis=0) + 0.0  # + 0.0: no "-0.0" for one species
        samples["taxonomy_db"] = tx["label"]
        samples["n_species"] = (spw > 0).sum(axis=0).reindex(order)
        samples["shannon"] = shannon.reindex(order)
        no_tax = [s for s in order if s not in tx["rel"].columns]
        samples.loc[no_tax, ["n_species", "shannon"]] = np.nan
        long = tx["long"]
        meta_join = samples[["batch"] + meta_cols]
        long = long.merge(meta_join, left_on="sample", right_index=True, how="left")
        lead = ["sample", "batch"] + meta_cols + ["taxon_id"] + tx["ranks"]
        long = long[lead + [c for c in long.columns if c not in lead]]
        long["_o"] = long["sample"].map({s: i for i, s in enumerate(order)})
        long = long.sort_values(["_o", "relative_abundance"], ascending=[True, False]).drop(columns="_o")
        tables["taxonomy_long"] = (long, False)
        for rank, wide in tx["wide"].items():
            tables[f"taxonomy_{rank.lower()}_wide"] = (wide[wide.sum(axis=1) > 0], True)
        tables["taxonomy_species_prevalence"] = (
            prevalence_table(spw, samples, args.group_by).rename_axis(sp), True)
        for s, sub in long.groupby("sample"):
            agg = sub.groupby(sp)["relative_abundance"].sum()
            species_by_sample[s] = {norm_taxon(n): (strip_gtdb(n), ab) for n, ab in agg.items() if n}
        summary["n_species"] = int((spw > 0).any(axis=1).sum())
        watch = load_watchlist(args.watchlist) if args.watchlist and Path(args.watchlist).exists() else None
        wl_wide = wl_long = None
        if watch is not None:
            wl_wide, wl_long = watchlist_table(tx, watch)
            tables["watchlist"] = (wl_wide, False)
        fig_taxonomy(figs, tx, [s for s in order if s in tx["rel"].columns],
                     [samples.at[s, args.group_by] if args.group_by else "" for s in order if s in tx["rel"].columns],
                     low, args.top_n, tx["label"])
        tcols = [s for s in order if s in tx["rel"].columns]
        fig_watchlist(figs, wl_wide, wl_long, tcols,
                      [samples.at[s, args.group_by] if args.group_by else "" for s in tcols], low)

    # ---- ARGs
    kframes = [b["kmarvel"] for b in batches if b["kmarvel"] is not None]
    ran = set().union(*[b["kmarvel_ran"] for b in batches])
    if kframes and ran:
        rules = load_rules(args.arg_priority)
        snps = load_snps(args.kmarvel_db)
        if snps is None:
            warn("no --kmarvel-db: tier 4 calls cannot be checked for resistance mutations")
        prevalence = load_prevalence(args.card_prevalence) if args.card_prevalence else None
        status = arg_status(samples, order, ran)
        samples["arg_status"] = [status[s] for s in order]
        not_run = [s for s in order if status[s] == NA_NOT_RUN]
        calls = pd.concat(kframes, ignore_index=True)
    if kframes and ran and calls.empty:
        warn("K-MARVEL found no ARG in any sample: ARG tables and figures are skipped")
        for col in [f"n_arg_calls_tier{t}" for t in (1, 2, 3, 4)] + ["acquired_genes", "acquired_rpkm"]:
            samples[col] = 0
        samples.loc[not_run, [f"n_arg_calls_tier{t}" for t in (1, 2, 3, 4)] + ["acquired_genes", "acquired_rpkm"]] = np.nan
        summary["n_tier1"], summary["n_tier1_genes"] = 0, 0
    elif kframes and ran:
        arg = annotate_args(calls, samples, rules, snps, prevalence, species_by_sample)
        mats = arg_tables(arg, samples, order, status, args.group_by, args.order_by)
        for t in (1, 2, 3, 4):
            samples[f"n_arg_calls_tier{t}"] = arg[arg["tier"] == t].groupby("sample").size().reindex(order).fillna(0).astype(int)
        acq = arg[(arg["tier"] <= 2) & arg["gene_representative"]]
        samples["acquired_genes"] = acq.groupby("sample")["gene"].nunique().reindex(order).fillna(0).astype(int)
        samples["acquired_rpkm"] = acq.groupby("sample")["rpkm"].sum().reindex(order).fillna(0)
        samples.loc[not_run, [f"n_arg_calls_tier{t}" for t in (1, 2, 3, 4)] + ["acquired_genes", "acquired_rpkm"]] = np.nan
        lead = ["sample", "batch"] + meta_cols + ["tier", "category", "gene", "aro_full_name", "aro_id", "gene_family",
                                                   "drug_class", "resistance_mechanism", "model", "identity",
                                                   "assembly_coverage", "supporting_reads", "median_kmer_depth",
                                                   "gene_length_nt", "rpm", "rpkm", "resistance_snp_status",
                                                   "resistance_snps_found", "kmer_blast_discordant",
                                                   "blast_predicted_aro_name", "gene_representative"]
        lead += [c for c in ("card_plasmid_associated", "candidate_hosts") if c in arg.columns]
        arg = arg.merge(samples[["batch"] + meta_cols], left_on="sample", right_index=True, how="left")
        arg["_o"] = arg["sample"].map({s: i for i, s in enumerate(order)})
        arg = arg.sort_values(["_o", "tier", "gene"]).drop(columns="_o")
        arg_out = arg[lead + ["tier_note", "mutation_data"]]
        tables["arg_long"] = (arg_out, False)
        for name, df in mats.items():
            tables[name] = (df, name != "group_arg_detection")
        hit1 = arg[arg["tier"] == 1]
        summary["n_tier1"] = int(hit1["sample"].nunique())
        summary["n_tier1_genes"] = int(hit1["gene"].nunique())
        fig_args(figs, arg[arg["gene_representative"]], mats, samples, order, groups, low, status, args.top_n, args.group_by, args.order_by)

    qc_lead = ["batch"] + meta_cols
    qc = samples[qc_lead + [c for c in samples.columns if c not in qc_lead]]
    tables = {"samples_qc": (qc.reset_index(), False), **tables}
    fig_qc(figs, samples, order, groups, low, args.min_clean_pairs)

    write_tables(outdir, tables)
    if not args.no_excel:
        write_excel(outdir, tables)
    versions = "\n".join(f"# {b['label']} ({b['root']})\n{b['versions']}" for b in batches if b["versions"])
    inputs = [f"{b['label']}: {b['root'].resolve()}" for b in batches]
    for flag, val in (("metadata", args.metadata), ("group-by", args.group_by), ("order-by", args.order_by),
                      ("kmarvel-db", args.kmarvel_db), ("card-prevalence", args.card_prevalence),
                      ("arg-priority", args.arg_priority), ("watchlist", args.watchlist),
                      ("min-clean-pairs", args.min_clean_pairs)):
        if val:
            inputs.append(f"--{flag} {val}")
    write_html(outdir, args.title, {
        "figures": figs, "tables": tables, "summary": summary, "versions": versions, "inputs": inputs,
        "min_pairs": args.min_clean_pairs, "meta_cols": meta_cols, "key_meta": key_meta,
        "group_by": args.group_by,
    })
    print(f"{len(order)} samples; {len(tables)} tables, {len(figs.records)} figures written to {outdir}")


if __name__ == "__main__":
    main()
