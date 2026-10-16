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
from summary_args import (METHOD_LABELS, annotate_args, annotate_rgi, arg_status, arg_tables,  # noqa: E402
                          compare_methods, load_kmarvel, load_prevalence, load_rgi, load_rules, load_snps)
from summary_figures import (FigureWriter, fig_args, fig_method_by_category, fig_method_upset,  # noqa: E402
                             fig_qc, fig_taxonomy, fig_watchlist, short_count)
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


def depth_classes(pairs, cuts):
    """Clean read pairs -> '<1M', '1M-10M', '10M-20M', '>=20M' style classes from the sorted cut points."""
    cuts = sorted(set(cuts))
    labels = [f"<{short_count(cuts[0])}"] + [f"{short_count(a)}-{short_count(b)}" for a, b in zip(cuts, cuts[1:])] \
        + [f">={short_count(cuts[-1])}"]
    idx = np.searchsorted(cuts, pairs.fillna(-1).values, side="right")
    return pd.Series([labels[i] if p >= 0 else "" for i, p in zip(idx, pairs.fillna(-1))], index=pairs.index)


def tier1_table(calls, comparison, primary):
    """Tier 1 calls of the primary method with the other method's evidence (when both ran)."""
    t1 = calls[calls["tier"] == 1].copy()
    keep = ["sample", "method", "gene", "aro_full_name", "gene_family", "category", "assembly_coverage",
            "supporting_reads", "rpkm"]
    keep += [c for c in ("rgi_kmer_hosts", "rgi_host_check", "candidate_hosts") if c in t1.columns]
    t1 = t1[keep]
    if comparison is not None:
        c = comparison["comp"]
        cols = ["sample", "gene", "detected_by"] + [x for x in ("rgi_gene_coverage", "rgi_kmer_hosts", "rgi_host_check")
                                                     if x in c.columns]
        t1 = t1.drop(columns=[x for x in ("rgi_kmer_hosts", "rgi_host_check") if x in t1.columns])
        t1 = t1.merge(c[cols], on=["sample", "gene"], how="left")
    return t1


def attention_flags(tier1, comparison, low_samples):
    """Short findings for the report header: [(kind, text, anchor)].

    Tier 1 genes come from the method comparison when several ARG methods ran (so a gene only one method found
    is listed too), else from the primary method's tier 1 calls.
    """
    flags = []
    rows = None
    if comparison is not None:
        rows = comparison["comp"]
        rows = rows[rows["tier"] == 1]
    elif tier1 is not None:
        rows = tier1.drop_duplicates(["sample", "gene"])
    if rows is not None and not rows.empty:
        for gene, sub in rows.groupby("gene", sort=True):
            n = len(sub)
            hosts = [h.split(" (")[0] for h in sub.get("rgi_kmer_hosts", pd.Series(dtype=str)).fillna("") if h]
            tops = sorted(((hosts.count(h), h) for h in set(hosts)), key=lambda x: (-x[0], x[1]))
            host_txt = ("; RGI k-mer top host: " + ", ".join(
                f"{h} ({c} sample{'s' if c > 1 else ''})" for c, h in tops[:3])) if tops else ""
            checks = set(sub.get("rgi_host_check", pd.Series(dtype=str)).fillna(""))
            text = f"{gene}: {n} sample{'s' if n > 1 else ''} ({', '.join(sub['sample'])}){host_txt}"
            if comparison is not None:
                by = sub["detected_by"].fillna("")
                n_all = int(by.str.contains(r"\+", regex=True).sum())
                single = sorted({b for b in by if "+" not in b})
            if "non_watchlist_host" in checks and "watchlist_host" not in checks:
                kind, note = "homolog", ("the RGI k-mer host is not a watchlist pathogen: possibly a related gene "
                                         "of a commensal")
            elif comparison is not None and n_all == 0:
                kind, note = "single", f"found by {' / '.join(single)}"
            elif comparison is not None:
                kind, note = "tier1", f"found by every method in {n_all} of {n} samples"
            else:
                kind, note = "tier1", ""
            flags.append((kind, text + (f" - {note}" if note else ""), "compare" if kind != "tier1" and
                          comparison is not None else "args"))
    if low_samples:
        flags.append(("low_depth", f"{len(low_samples)} low-depth sample{'s' if len(low_samples) > 1 else ''}: "
                                   f"{', '.join(low_samples)} (absence of a gene or taxon there is not informative)",
                      "qc"))
    return flags


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
    p.add_argument("--rgi-min-coverage", type=float, default=80,
                   help="RGI bwt: a gene is detected when reads cover at least this %% of its best-covered "
                        "reference sequence (default 80). RGI itself reports every gene with any mapped read")
    p.add_argument("--arg-primary", choices=("auto", "kmarvel", "rgi"), default="auto",
                   help="read-level ARG method behind the ARG matrices and figures (default auto: K-MARVEL when "
                        "present, else RGI); every method's calls are in arg_long.tsv")
    p.add_argument("--depth-marks", default="10000000,20000000",
                   help="comma-separated clean read-pair depths marked in the QC figure and counted in the report "
                        "header (default 10000000,20000000; '' for none)")
    p.add_argument("--taxa-per-sample", type=int, default=5,
                   help="the taxonomy heatmaps show the union of each sample's N most abundant taxa (default 5)")
    p.add_argument("--top-n", type=int, default=40, help="maximum rows in the ARG heatmap (default 40)")
    p.add_argument("--title", default="BugBuster results summary", help="HTML report title")
    p.add_argument("--no-excel", action="store_true", help="do not write the Excel workbook")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.order_by and not args.group_by:
        warn("--order-by without --group-by orders all samples as one series")
    try:
        depth_marks = sorted({int(float(x)) for x in str(args.depth_marks).split(",") if x.strip()})
    except ValueError:
        die(f"--depth-marks: '{args.depth_marks}' is not a comma-separated list of read-pair counts")

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
        rgi, rgi_kmer, rgi_ran, rgi_kmer_ran = load_rgi(root)
        names = set()
        if reads is not None:
            names |= set(reads["sample"])
        if taxonomy is not None:
            names |= set(taxonomy["otu"].columns)
        names |= kran | rgi_ran
        if not names:
            die(f"{root}: no BugBuster outputs found (expected 01_quality_control/, 02_taxonomy/ or 05_arg_prediction/)")
        for s in names:
            if s in seen:
                die(f"sample {s} is in both {seen[s]} and {root}; sample names must be unique across results directories")
            seen[s] = str(root)
        batches.append({"label": label, "root": root, "reads": reads, "fastp": fastp, "taxonomy": taxonomy,
                        "kmarvel": kmarvel, "kmarvel_ran": kran, "rgi": rgi, "rgi_kmer": rgi_kmer, "rgi_ran": rgi_ran,
                        "rgi_kmer_ran": rgi_kmer_ran, "samples": names, "versions": load_versions(root)})
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
        samples["depth_class"] = depth_classes(samples["clean_pairs"], [args.min_clean_pairs] + depth_marks)
    else:
        warn("no Reads_report.csv: depth flags and ARG normalisation are unavailable")
    lengths = pd.concat([b["fastp"] for b in batches])
    if not lengths.empty:
        samples = samples.join(lengths.set_index("sample"), how="left")
    low = set(samples.index[samples["low_depth"].fillna(False).astype(bool)]) if "low_depth" in samples else set()

    tables = {}
    figs = FigureWriter(outdir)
    summary = {"n_samples": len(order), "n_batches": len(batches), "n_low_depth": len(low),
               "low_depth_samples": [s for s in order if s in low]}
    if "clean_pairs" in samples:
        summary["depth_marks"] = [(m, int((samples["clean_pairs"] >= m).sum())) for m in depth_marks]

    watch = load_watchlist(args.watchlist) if args.watchlist and Path(args.watchlist).exists() else None

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
        wl_wide = wl_long = None
        if watch is not None:
            wl_wide, wl_long = watchlist_table(tx, watch)
            tables["watchlist"] = (wl_wide, False)
        tcols = [s for s in order if s in tx["rel"].columns]
        tgroups = [samples.at[s, args.group_by] if args.group_by else "" for s in tcols]
        fig_taxonomy(figs, tx, tcols, tgroups, low, args.taxa_per_sample, tx["label"])
        fig_watchlist(figs, wl_wide, wl_long, tcols, tgroups, low, Path(args.watchlist).name)

    # ---- ARGs (read level): K-MARVEL calls, and RGI bwt genes at or above --rgi-min-coverage
    methods_ran = {}
    kran = set().union(*[b["kmarvel_ran"] for b in batches])
    rgi_ran = set().union(*[b["rgi_ran"] for b in batches])
    if kran:
        methods_ran["kmarvel"] = kran
    if rgi_ran:
        methods_ran["rgi"] = rgi_ran
    comparison, rgi_all, calls_by = None, None, {}
    if methods_ran:
        primary = args.arg_primary
        if primary == "auto":
            primary = "kmarvel" if "kmarvel" in methods_ran else "rgi"
        elif primary not in methods_ran:
            die(f"--arg-primary {primary}: no {METHOD_LABELS[primary]} output in the results directories")
        rules = load_rules(args.arg_priority)
        prevalence = load_prevalence(args.card_prevalence) if args.card_prevalence else None
        if "kmarvel" in methods_ran:
            snps = load_snps(args.kmarvel_db)
            if snps is None:
                warn("no --kmarvel-db: tier 4 calls cannot be checked for resistance mutations")
            kcalls = pd.concat([b["kmarvel"] for b in batches if b["kmarvel"] is not None], ignore_index=True)
            calls_by["kmarvel"] = annotate_args(kcalls, samples, rules, snps, prevalence, species_by_sample)
        if "rgi" in methods_ran:
            alleles = pd.concat([b["rgi"] for b in batches if b["rgi"] is not None], ignore_index=True)
            kmer = pd.concat([b["rgi_kmer"] for b in batches if b["rgi_kmer"] is not None], ignore_index=True)
            kmer_ran = set().union(*[b["rgi_kmer_ran"] for b in batches])
            rgi_all = annotate_rgi(alleles, kmer, kmer_ran, samples, rules, prevalence, species_by_sample, watch,
                                   args.rgi_min_coverage)
            calls_by["rgi"] = rgi_all[rgi_all["rgi_detected"]].copy() if not rgi_all.empty else rgi_all
        summary["arg_primary"] = METHOD_LABELS[primary]
        summary["arg_methods"] = [METHOD_LABELS[m] for m in methods_ran]
        status = arg_status(samples, order, methods_ran[primary])
        samples["arg_status"] = [status[s] for s in order]
        not_run = [s for s in order if status[s] == NA_NOT_RUN]
        count_cols = [f"n_arg_calls_tier{t}" for t in (1, 2, 3, 4)] + ["acquired_genes", "acquired_rpkm"]
        arg = calls_by[primary]
        if arg.empty:
            warn(f"{METHOD_LABELS[primary]} found no ARG in any sample: ARG tables and figures are skipped")
            for col in count_cols:
                samples[col] = 0
            samples.loc[not_run, count_cols] = np.nan
            summary["n_tier1"], summary["n_tier1_genes"] = 0, 0
        else:
            mats = arg_tables(arg, samples, order, status, args.group_by, args.order_by)
            for t in (1, 2, 3, 4):
                samples[f"n_arg_calls_tier{t}"] = arg[arg["tier"] == t].groupby("sample").size().reindex(order).fillna(0).astype(int)
            acq = arg[(arg["tier"] <= 2) & arg["gene_representative"]]
            samples["acquired_genes"] = acq.groupby("sample")["gene"].nunique().reindex(order).fillna(0).astype(int)
            samples["acquired_rpkm"] = acq.groupby("sample")["gene_rpkm"].sum().reindex(order).fillna(0)
            samples.loc[not_run, count_cols] = np.nan
            for name, df in mats.items():
                tables[name] = (df, name != "group_arg_detection")
            hit1 = arg[arg["tier"] == 1]
            summary["n_tier1"] = int(hit1["sample"].nunique())
            summary["n_tier1_genes"] = int(hit1["gene"].nunique())
            fig_args(figs, arg[arg["gene_representative"]], mats, samples, order, groups, low, status, args.top_n,
                     METHOD_LABELS[primary])
        # arg_long: every method's calls, the primary method first
        lead = ["sample", "batch"] + meta_cols + ["method", "tier", "category", "gene", "aro_full_name", "aro_id",
                                                   "gene_family", "drug_class", "resistance_mechanism", "model",
                                                   "identity", "assembly_coverage", "supporting_reads",
                                                   "median_kmer_depth", "gene_length_nt", "rpm", "rpkm", "gene_rpkm",
                                                   "resistance_snp_status", "resistance_snps_found",
                                                   "kmer_blast_discordant", "blast_predicted_aro_name",
                                                   "gene_representative"]
        frames = [calls_by[m] for m in [primary] + [m for m in calls_by if m != primary] if not calls_by[m].empty]
        if frames:
            allc = pd.concat(frames, ignore_index=True)
            extra = [c for c in ("card_plasmid_associated", "candidate_hosts", "rgi_mapped_reads", "rgi_best_ref",
                                 "rgi_gene_coverage", "rgi_kmer_hosts", "rgi_host_check") if c in allc.columns]
            allc = allc.merge(samples[["batch"] + meta_cols], left_on="sample", right_index=True, how="left")
            allc["_o"] = allc["sample"].map({s: i for i, s in enumerate(order)})
            allc["_m"] = (allc["method"] != primary).astype(int)
            allc = allc.sort_values(["_o", "_m", "tier", "gene"], kind="stable").drop(columns=["_o", "_m"])
            tables["arg_long"] = (allc[lead + extra + ["tier_note", "mutation_data"]], False)
        if len(calls_by) > 1:
            comp, overlap, gene_sets = compare_methods(calls_by, rgi_all, methods_ran)
            if not comp.empty:
                lead_c = ["sample", "batch"] + meta_cols
                comp = comp.merge(samples[["batch"] + meta_cols], left_on="sample", right_index=True, how="left")
                comp["_o"] = comp["sample"].map({s: i for i, s in enumerate(order)})
                comp = comp.sort_values(["_o", "tier", "gene"], kind="stable").drop(columns="_o")
                comp = comp[lead_c + [c for c in comp.columns if c not in lead_c]]
                tables["arg_method_comparison"] = (comp, False)
                tables["arg_method_overlap"] = (overlap, False)
                common = set.intersection(*[set(v) for v in methods_ran.values()])
                for m in calls_by:
                    n = comp[comp[f"detected_{m}"]].groupby("sample")["gene"].nunique()
                    samples[f"n_genes_{m}"] = n.reindex(order).fillna(0).astype(int)
                both = comp[comp[[f"detected_{m}" for m in calls_by]].all(axis=1)]
                samples["n_genes_all_methods"] = both.groupby("sample")["gene"].nunique().reindex(order).fillna(0).astype(int)
                outside = [s for s in order if s not in common]
                samples.loc[outside, [f"n_genes_{m}" for m in calls_by] + ["n_genes_all_methods"]] = np.nan
                summary["n_tier1_genes_all_methods"] = int(both.loc[both["tier"] == 1, "gene"].nunique())
                comparison = {"comp": comp, "overlap": overlap, "gene_sets": gene_sets,
                              "methods": [m for m in ("kmarvel", "rgi") if m in calls_by]}
                fig_method_upset(figs, comp, gene_sets, comparison["methods"])
                fig_method_by_category(figs, comp, comparison["methods"])

    prim_arg = None
    if methods_ran:
        prim_arg = calls_by[primary]
    tier1_view = None
    if prim_arg is not None and not prim_arg.empty:
        tier1_view = tier1_table(calls_by[primary], comparison, primary)
    summary["attention"] = attention_flags(tier1_view, comparison, summary["low_depth_samples"])

    qc_lead = ["batch"] + meta_cols
    qc = samples[qc_lead + [c for c in samples.columns if c not in qc_lead]]
    tables = {"samples_qc": (qc.reset_index(), False), **tables}
    fig_qc(figs, samples, order, groups, low, args.min_clean_pairs, depth_marks)

    write_tables(outdir, tables)
    if not args.no_excel:
        write_excel(outdir, tables)
    versions = "\n".join(f"# {b['label']} ({b['root']})\n{b['versions']}" for b in batches if b["versions"])
    inputs = [f"{b['label']}: {b['root'].resolve()}" for b in batches]
    for flag, val in (("metadata", args.metadata), ("group-by", args.group_by), ("order-by", args.order_by),
                      ("kmarvel-db", args.kmarvel_db), ("card-prevalence", args.card_prevalence),
                      ("arg-priority", args.arg_priority), ("watchlist", args.watchlist),
                      ("min-clean-pairs", args.min_clean_pairs), ("rgi-min-coverage", args.rgi_min_coverage),
                      ("arg-primary", args.arg_primary), ("depth-marks", args.depth_marks),
                      ("taxa-per-sample", args.taxa_per_sample)):
        if val:
            inputs.append(f"--{flag} {val}")
    write_html(outdir, args.title, {
        "figures": figs, "tables": tables, "summary": summary, "versions": versions, "inputs": inputs,
        "min_pairs": args.min_clean_pairs, "meta_cols": meta_cols, "key_meta": key_meta,
        "group_by": args.group_by, "comparison": comparison, "tier1_view": tier1_view,
        "rgi_min_coverage": args.rgi_min_coverage, "watchlist": watch, "watchlist_path": args.watchlist,
        "arg_priority": args.arg_priority,
    })
    print(f"{len(order)} samples; {len(tables)} tables, {len(figs.records)} figures written to {outdir}")


if __name__ == "__main__":
    main()
