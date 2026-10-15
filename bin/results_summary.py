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
import html
import json
import math
import re
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_ASSETS = SCRIPT_DIR.parent / "assets"

TIER_LABELS = {
    1: "Tier 1: clinical priority",
    2: "Tier 2: other acquired",
    3: "Tier 3: intrinsic / efflux / regulator",
    4: "Tier 4: mutation dependent",
}
NA_LOW_DEPTH = "NA_low_depth"
NA_NOT_RUN = "NA_not_run"
GTDB_SUFFIX = re.compile(r"_[A-Z]+$")
ACCESSION = re.compile(r"(GC[AF]_\d+\.\d+)")
KMARVEL_COLUMNS = [
    "aro_id", "aro_full_name", "gene_family", "drug_class", "resistance_mechanism",
    "model", "identity", "assembly_coverage", "supporting_reads", "median_kmer_depth",
    "blast_predicted_aro_name", "blast_predicted_aro_id", "mutation_data",
    "reference_sequence",
]


# --------------------------------------------------------------------------- helpers

def die(msg):
    sys.exit(f"ERROR: {msg}")


def warn(msg):
    print(f"WARNING: {msg}", file=sys.stderr)


def natural_key(value):
    """Sort key that orders V2 before V10."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(value))]


def strip_gtdb(name):
    """Remove GTDB placeholder suffixes for display and matching: Enterococcus_B faecium -> Enterococcus faecium."""
    if not isinstance(name, str):
        return name
    return " ".join(GTDB_SUFFIX.sub("", word) for word in name.split(" "))


def norm_taxon(name):
    return strip_gtdb(name).strip().lower() if isinstance(name, str) else ""


def batch_label(path):
    p = Path(path).resolve()
    return p.parent.name if p.name == "results" else p.name


# --------------------------------------------------------------------------- loading

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


def load_kmarvel(root):
    kdir = root / "05_arg_prediction" / "read_level" / "kmarvel"
    if not kdir.is_dir():
        return None, set()
    frames, ran = [], set()
    for path in sorted(kdir.glob("*/*_kmarvel.tsv")):
        sample = path.name[: -len("_kmarvel.tsv")]
        ran.add(sample)
        header = pd.read_csv(path, sep="\t", nrows=0).columns
        missing = [c for c in KMARVEL_COLUMNS if c not in header]
        if missing:
            die(f"{path} lacks columns: {', '.join(missing)}")
        # read_data and the assembled sequences grow with depth: read only what is used
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, usecols=KMARVEL_COLUMNS)
        if df.empty:
            continue
        df.insert(0, "sample", sample)
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["sample"] + KMARVEL_COLUMNS), ran
    return pd.concat(frames, ignore_index=True), ran


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


def load_snps(kmarvel_db):
    """CARD snps.txt -> {aro_number: [set of required mutations, ...]} (single, multiple, co-dependent variants)."""
    if not kmarvel_db:
        return None
    path = Path(kmarvel_db) / "snps.txt"
    if not path.exists():
        warn(f"{path} not found: tier 4 calls get resistance_snp_status 'not_checked'")
        return None
    snps = pd.read_csv(path, sep="\t", dtype=str)
    rules = {}
    for acc, ptype, muts in snps[["Accession", "Parameter Type", "Mutations"]].itertuples(index=False):
        key = str(acc).strip()
        rules.setdefault(key, [])
        if "variant" not in str(ptype) and "nonsense" not in str(ptype):
            continue  # frameshifts cannot be checked from K-MARVEL's substitution list
        required = {normalise_mutation(m) for m in str(muts).split(",") if m.strip()}
        if required:
            rules[key].append(required)
    return rules


def normalise_mutation(token):
    return token.strip().replace("Ter", "*")


def load_prevalence(path):
    """WildCARD card_prevalence.txt -> {ARO: (set of normalised host species, plasmid-associated bool)}."""
    if not path:
        return None
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = {"ARO Accession", "Pathogen", "NCBI Plasmid"}
    if not need <= set(df.columns):
        die(f"{path} is not a CARD prevalence table (needs {sorted(need)})")
    df["plasmid"] = pd.to_numeric(df["NCBI Plasmid"], errors="coerce").fillna(0) > 0
    out = {}
    for aro, grp in df.groupby("ARO Accession"):
        out[aro] = ({norm_taxon(p) for p in grp["Pathogen"]}, bool(grp["plasmid"].any()))
    return out


def load_rules(path):
    rules = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    need = {"match", "pattern", "tier", "category"}
    if not need <= set(rules.columns):
        die(f"{path} needs columns {sorted(need)}")
    if "note" not in rules.columns:
        rules["note"] = ""
    compiled = []
    for r in rules.itertuples(index=False):
        if r.match not in ("name", "family"):
            die(f"{path}: match must be 'name' or 'family', got '{r.match}'")
        compiled.append((r.match, re.compile(r.pattern), int(r.tier), r.category, r.note))
    return compiled


def load_watchlist(path):
    wl = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    if not {"taxon", "rank"} <= set(wl.columns):
        die(f"{path} needs columns taxon, rank")
    if "note" not in wl.columns:
        wl["note"] = ""
    return wl


# --------------------------------------------------------------------------- ARG annotation

def gene_name(aro, family):
    """Collapse allele numbers, which are unstable in read-level calls (KPC-2 -> KPC)."""
    if aro.startswith("OXA-"):
        m = re.search(r"(OXA-\d+-like)", family)
        return m.group(1) if m else aro
    m = re.match(r"^MCR-(\d+)", aro)
    if m:
        return f"MCR-{m.group(1)}"
    m = re.match(r"^(AAC\(6'\)-Ib-cr)\d*$", aro)
    if m:
        return m.group(1)
    m = re.match(r"^(Cfx[A-Z]|Qnr[A-Z]{1,2})\d+$", aro)
    if m:
        return m.group(1)
    m = re.match(r"^([A-Za-z]{2,6}(?:-[A-Z])?)-\d+[a-z]?(?:\.\d+)?$", aro)
    if m:
        return m.group(1)
    return aro


def snp_status(aro_id, mutation_data, snps):
    if snps is None:
        return "not_checked", ""
    key = aro_id.replace("ARO:", "")
    defs = snps.get(key)
    if not defs:
        return "no_snp_definition", ""
    seen = {normalise_mutation(m) for m in re.split(r"[;,]", mutation_data or "") if m.strip()}
    found = [",".join(sorted(req)) for req in defs if req <= seen]
    if found:
        return "card_snp_present", ";".join(sorted(set(found)))
    return "no_card_snp", ""


def assign_tier(row, rules, snps):
    name, family = row["aro_full_name"], row["gene_family"]
    key = row["aro_id"].replace("ARO:", "")
    mutation_model = (
        row["model"] == "protein variant model"
        or (snps is not None and bool(snps.get(key)))
        or re.search(r"mutant|mutation", name, re.I) is not None
    )
    if mutation_model:
        return 4, "mutation_dependent", "resistance needs a known mutation (see resistance_snp_status)"
    for match, pattern, tier, category, note in rules:
        target = name if match == "name" else family
        if pattern.search(target):
            return tier, category, note
    mech = row["resistance_mechanism"].lower()
    if "efflux" in mech or "permeability" in mech or row["model"] == "protein overexpression model":
        return 3, "efflux_permeability", "efflux pump, porin or regulator (often chromosomal)"
    return 2, "acquired_unreviewed", "no tier rule matched; review and add a rule to the tier table"


def annotate_args(arg, samples, rules, snps, prevalence, species_by_sample):
    if arg.empty:
        return arg
    arg = arg.copy()
    for col in ("identity", "assembly_coverage", "supporting_reads", "median_kmer_depth"):
        arg[col] = pd.to_numeric(arg[col], errors="coerce")
    arg["gene_length_nt"] = arg["reference_sequence"].str.replace(r"[^A-Za-z]", "", regex=True).str.len() * 3
    arg.loc[arg["gene_length_nt"] == 0, "gene_length_nt"] = np.nan
    clean = arg["sample"].map(samples["clean_pairs"]) if "clean_pairs" in samples else np.nan
    arg["rpm"] = arg["supporting_reads"] / clean * 1e6
    arg["rpkm"] = arg["rpm"] / (arg["gene_length_nt"] / 1000)
    tiers = arg.apply(lambda r: assign_tier(r, rules, snps), axis=1, result_type="expand")
    arg["tier"], arg["category"], arg["tier_note"] = tiers[0].astype(int), tiers[1], tiers[2]
    arg["gene"] = [gene_name(a, f) for a, f in zip(arg["aro_full_name"], arg["gene_family"])]
    status = [
        snp_status(a, m, snps) if t == 4 else ("", "")
        for a, m, t in zip(arg["aro_id"], arg["mutation_data"], arg["tier"])
    ]
    arg["resistance_snp_status"] = [s[0] for s in status]
    arg["resistance_snps_found"] = [s[1] for s in status]
    blast = arg["blast_predicted_aro_id"].fillna("")
    arg["kmer_blast_discordant"] = (blast != "") & (blast != arg["aro_id"])
    # alleles of one gene in one sample mostly share their reads (TEM-1 and TEM-256 with 16032 and 16029
    # pairs): the gene is represented by its best-supported allele, never by the sum
    arg = arg.sort_values(["sample", "gene", "supporting_reads"], ascending=[True, True, False])
    arg["gene_representative"] = ~arg.duplicated(["sample", "gene"])
    if prevalence is not None:
        hosts_col, plasmid_col = [], []
        for s, aro in zip(arg["sample"], arg["aro_id"]):
            hosts, plasmid = prevalence.get(aro, (set(), None))
            plasmid_col.append("" if plasmid is None else plasmid)
            detected = species_by_sample.get(s, {})
            hits = sorted(
                ((ab, disp) for norm, (disp, ab) in detected.items() if norm in hosts),
                reverse=True,
            )
            hosts_col.append("; ".join(f"{d} ({ab:.2g}%)" for ab, d in hits[:5]))
        arg["card_plasmid_associated"] = plasmid_col
        arg["candidate_hosts"] = hosts_col
    return arg


# --------------------------------------------------------------------------- summary tables

def sample_order(samples, group_by, order_by):
    keys = []
    for s in samples.index:
        g = samples.at[s, group_by] if group_by else ""
        o = samples.at[s, order_by] if order_by else ""
        keys.append((natural_key(g), natural_key(o), natural_key(s), s))
    return [k[-1] for k in sorted(keys)]


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


def arg_status(samples, order, ran):
    """Per sample: ok, NA_low_depth (absence not informative) or NA_not_run (no K-MARVEL output)."""
    low = samples["low_depth"].fillna(False).astype(bool) if "low_depth" in samples else pd.Series(False, samples.index)
    return {s: NA_NOT_RUN if s not in ran else (NA_LOW_DEPTH if low[s] else "ok") for s in order}


def arg_tables(arg, samples, order, status, group_by, order_by):
    """Wide ARG matrices; a cell for a gene not found is 0, NA_low_depth or NA_not_run."""
    out = {}
    arg = arg[arg["gene_representative"]]
    cols = [s for s in order]

    def fill_missing(mat):
        mat = mat.reindex(columns=cols)
        mat = mat.astype(object)
        for s in cols:
            fill = 0 if status[s] == "ok" else status[s]
            mat[s] = mat[s].where(mat[s].notna(), fill)
        return mat

    gene_meta = arg.groupby("gene").agg(tier=("tier", "min"), category=("category", "first"),
                                        gene_family=("gene_family", "first"))
    rpkm = arg.pivot_table(index="gene", columns="sample", values="rpkm", aggfunc="sum")
    rpkm_numeric = rpkm.reindex(columns=cols)
    gene_rpkm = fill_missing(rpkm)
    gene_rpkm.insert(0, "category", gene_meta.loc[gene_rpkm.index, "category"])
    gene_rpkm.insert(0, "tier", gene_meta.loc[gene_rpkm.index, "tier"])
    gene_rpkm = gene_rpkm.sort_values(["tier"], kind="stable")
    out["arg_gene_rpkm"] = gene_rpkm

    fam = arg.pivot_table(index="gene_family", columns="sample", values="rpkm", aggfunc="sum")
    out["arg_family_rpkm"] = fill_missing(fam)

    pres = (rpkm_numeric.reindex(columns=cols).notna() & (rpkm_numeric.reindex(columns=cols) > 0)).astype(int)
    pres = pres.astype(object)
    for s in cols:
        if status[s] != "ok":
            pres[s] = pres[s].where(pres[s] == 1, status[s])
    pres.insert(0, "tier", gene_meta.loc[pres.index, "tier"])
    out["arg_presence"] = pres.sort_values("tier", kind="stable")

    counted = arg[(arg["tier"] < 4) | (arg["resistance_snp_status"] == "card_snp_present")].copy()
    counted["tier_group"] = np.where(counted["tier"] <= 2, "acquired (tiers 1-2)",
                                     np.where(counted["tier"] == 3, "intrinsic/efflux (tier 3)",
                                              "CARD mutation present (tier 4)"))
    counted["drug_class"] = counted["drug_class"].str.split(";")
    counted = counted.explode("drug_class")
    counted["drug_class"] = counted["drug_class"].str.strip()
    counted = counted[counted["drug_class"] != ""]
    dc = counted.pivot_table(index=["tier_group", "drug_class"], columns="sample", values="rpkm", aggfunc="sum")
    out["arg_drug_class_rpkm"] = fill_missing(dc) if not dc.empty else dc

    t1 = rpkm_numeric.loc[gene_meta.loc[rpkm_numeric.index, "tier"] == 1]
    panel = fill_missing(t1)
    panel.insert(0, "category", gene_meta.loc[panel.index, "category"])
    out["arg_priority_panel"] = panel

    if group_by:
        rows = []
        sel = arg[arg["tier"] <= 2]
        for g, members in samples.loc[cols].groupby(group_by, sort=False):
            mem = [s for s in cols if s in members.index]
            labels = {s: (samples.at[s, order_by] if order_by else s) for s in mem}
            detected_genes = sel[sel["sample"].isin(mem)]["gene"].unique()
            for gene in sorted(detected_genes, key=lambda x: (gene_meta.at[x, "tier"], x.lower())):
                hit = set(sel[(sel["gene"] == gene) & sel["sample"].isin(mem)]["sample"])
                pattern, evaluable = [], 0
                for s in mem:
                    if s in hit:
                        mark = "+"
                    elif status[s] != "ok":
                        mark = "NA"
                    else:
                        mark = "-"
                    evaluable += status[s] == "ok" or s in hit
                    pattern.append(f"{labels[s]}:{mark}")
                det = [labels[s] for s in mem if s in hit]
                rows.append({
                    group_by: g, "gene": gene, "tier": gene_meta.at[gene, "tier"],
                    "category": gene_meta.at[gene, "category"],
                    "n_detected": len(hit), "n_evaluable": evaluable, "n_samples": len(mem),
                    "first_detected": det[0] if det else "", "last_detected": det[-1] if det else "",
                    "pattern": " ".join(pattern),
                })
        out["group_arg_detection"] = pd.DataFrame(rows)
    return out


# --------------------------------------------------------------------------- figures

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


# --------------------------------------------------------------------------- outputs

COLUMN_DOCS = {
    "samples_qc": [
        ("sample / batch", "sample name and the results directory it came from"),
        ("raw_pairs, qc_pass_pairs, clean_pairs", "read pairs before QC, after fastp, after host/contaminant removal "
                                                  "(Reads_report.csv counts single reads; these are reads / 2)"),
        ("pct_removed_qc, pct_removed_decontam, pct_clean", "share of raw read pairs removed by fastp, removed by "
                                                            "host/contaminant filtering, and kept"),
        ("low_depth", "clean_pairs below --min-clean-pairs; absence calls in this sample are not informative"),
        ("classified_fraction", "fraction of reads the taxonomic profiler classified (sylph: skewed upward when "
                                "--sylph_read_seq_id is not set and many reads come from genomes absent from the database)"),
        ("n_species, shannon", "species detected and Shannon index (natural log) of species relative abundance"),
        ("arg_status", "ok, NA_low_depth or NA_not_run (K-MARVEL output missing)"),
        ("n_arg_calls_tier1..4", "K-MARVEL calls per tier (tier 4 = mutation dependent; only rows with "
                                 "resistance_snp_status card_snp_present carry a known resistance mutation)"),
        ("acquired_genes, acquired_rpkm", "distinct tier 1-2 genes and their summed RPKM"),
    ],
    "taxonomy_long": [
        ("relative_abundance", "% of the profiled (classified) community; sylph: taxonomic abundance"),
        ("est_read_pairs", "read pairs assigned to the taxon (sylph: sequence_abundance / 100 x clean read pairs)"),
        ("sequence_abundance", "sylph: % of all reads, unclassified included"),
        ("adjusted_ani, ani_5_95, true_cov, eff_lambda", "sylph containment ANI, its 5-95% interval (low coverage), "
                                                         "estimated coverage and effective lambda"),
        ("species_call_uncertain", "the lower bound of sylph's ANI interval is below 95%: the genome present may be "
                                   "a different species of the same genus"),
    ],
    "arg_long": [
        ("tier, category, tier_note", "reporting tier from the tier table (1 clinical priority, 2 other acquired, "
                                      "3 intrinsic/efflux/regulator, 4 mutation dependent)"),
        ("gene", "allele numbers collapsed (KPC-2 -> KPC; OXA alleles -> their CARD OXA subgroup)"),
        ("supporting_reads", "K-MARVEL supporting read pairs"),
        ("gene_length_nt", "CARD reference protein length x 3"),
        ("rpm, rpkm", "read pairs per million clean read pairs; RPM per kb of reference gene"),
        ("resistance_snp_status", "tier 4: card_snp_present (a CARD snps.txt resistance mutation is among "
                                  "K-MARVEL's mutation_data), no_card_snp, no_snp_definition, not_checked (no "
                                  "--kmarvel-db). CARD lists some common polymorphisms (E. coli MarR Y137H), so "
                                  "card_snp_present is necessary, not sufficient, evidence of resistance"),
        ("gene_representative", "True for the best-supported allele of each gene in a sample; matrices and "
                                "figures use these rows only, because alleles of one gene share their reads"),
        ("kmer_blast_discordant", "K-MARVEL's k-mer assignment and its BLAST best hit name different AROs"),
        ("card_plasmid_associated", "CARD prevalence lists this ARO on NCBI plasmids"),
        ("candidate_hosts", "species detected in the same sample that CARD prevalence lists as carriers of this "
                            "ARO, with their abundance. Co-occurrence only: reads are not linked to a host"),
    ],
    "arg matrices": [
        ("arg_gene_rpkm / arg_family_rpkm", "RPKM of the gene's best-supported allele / summed over a CARD family's "
                                            "genes; 0 = not detected in a sample of "
                                            "adequate depth; NA_low_depth / NA_not_run = unknown"),
        ("arg_presence", "1 detected, 0 not detected, NA_low_depth / NA_not_run unknown"),
        ("arg_drug_class_rpkm", "RPKM per drug class, split into acquired (tiers 1-2), intrinsic/efflux (tier 3) "
                                "and tier 4 calls with a CARD resistance mutation; a multi-class gene counts in "
                                "each class"),
        ("arg_priority_panel", "tier 1 genes only"),
        ("group_arg_detection", "per group and tier 1-2 gene: detections, evaluable samples and the detection "
                                "pattern in --order-by order (+ detected, - not detected, NA unknown)"),
    ],
}


def write_tables(outdir, tables):
    tdir = outdir / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    for name, (df, index) in tables.items():
        df.to_csv(tdir / f"{name}.tsv", sep="\t", index=index)


def write_excel(outdir, tables):
    engine = None
    for mod in ("openpyxl", "xlsxwriter"):
        try:
            __import__(mod)
            engine = mod
            break
        except ImportError:
            continue
    if engine is None:
        warn("neither openpyxl nor xlsxwriter is installed: no Excel workbook written")
        return None
    path = outdir / "results_summary.xlsx"
    readme = []
    for section, cols in COLUMN_DOCS.items():
        for col, desc in cols:
            readme.append({"table": section, "column": col, "description": desc})
    with pd.ExcelWriter(path, engine=engine) as xw:
        pd.DataFrame(readme).to_excel(xw, sheet_name="README", index=False)
        for name, (df, index) in tables.items():
            if len(df) > 1_000_000:
                warn(f"{name} has {len(df)} rows: too many for Excel, see tables/{name}.tsv")
                continue
            out = df.reset_index() if index else df
            out.to_excel(xw, sheet_name=name[:31], index=False, freeze_panes=(1, 0))
    return path


def df_to_html(df, max_rows=500, index=False, table_id=""):
    if df is None or df.empty:
        return "<p class='muted'>No rows.</p>"
    d = df.reset_index() if index else df
    note = ""
    if len(d) > max_rows:
        note = f"<p class='muted'>First {max_rows} of {len(d)} rows; the full table is in the TSV/Excel output.</p>"
        d = d.head(max_rows)

    def fmt(v):
        if isinstance(v, (float, np.floating)):
            if np.isnan(v):
                return ""
            if v == 0:
                return "0"
            return f"{v:.3g}" if abs(v) < 1000 else f"{v:,.0f}"
        if isinstance(v, (int, np.integer)) and not isinstance(v, bool):
            return f"{v:,}" if abs(v) >= 10000 else str(v)
        return html.escape(str(v))

    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in d.columns)
    body = []
    for row in d.itertuples(index=False):
        cells = []
        for v in row:
            raw = v if isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool) else ""
            sort = f" data-v='{raw}'" if raw != "" and not (isinstance(raw, float) and np.isnan(raw)) else ""
            cls = " class='na'" if isinstance(v, str) and v.startswith("NA_") else ""
            cells.append(f"<td{sort}{cls}>{fmt(v)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return (f"{note}<div class='tablewrap'><input class='filter' placeholder='Filter rows…' "
            f"data-target='{table_id}'><table class='sortable' id='{table_id}'><thead><tr>{head}</tr></thead>"
            f"<tbody>{''.join(body)}</tbody></table></div>")


def svg_inline(path):
    text = path.read_text()
    text = text[text.find("<svg"):]
    text = re.sub(r'<svg([^>]*?) width="[^"]+" height="[^"]+"', r'<svg\1 class="fig"', text, count=1)
    return text


HTML_CSS = """
:root{--bg:#ffffff;--fg:#1d2329;--muted:#5f6b76;--line:#d9dee3;--card:#f6f8fa;--accent:#2a7ab9;--warn:#c0392b;--figbg:#ffffff}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#15191d;--fg:#e4e8ec;--muted:#9aa5b1;--line:#2c343c;--card:#1d2328;--accent:#6fb1e6;--warn:#ef7b6c;--figbg:#f4f4f4}}
:root[data-theme="dark"]{--bg:#15191d;--fg:#e4e8ec;--muted:#9aa5b1;--line:#2c343c;--card:#1d2328;--accent:#6fb1e6;--warn:#ef7b6c;--figbg:#f4f4f4}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
header{padding:28px 16px 12px;max-width:1200px;margin:0 auto}
main{max-width:1200px;margin:0 auto;padding:0 16px 48px}
nav{position:sticky;top:0;background:var(--bg);border-bottom:1px solid var(--line);padding:8px 16px;z-index:2;overflow-x:auto;white-space:nowrap}
nav a{color:var(--accent);margin-right:16px;text-decoration:none;font-size:14px}
h1{margin:0 0 6px;font-size:26px}h2{margin:40px 0 10px;font-size:20px;border-bottom:1px solid var(--line);padding-bottom:6px}
h3{margin:24px 0 8px;font-size:16px}
.muted{color:var(--muted);font-size:13px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:10px;margin:16px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px 12px}
.card b{display:block;font-size:22px}.card span{color:var(--muted);font-size:13px}
.card.warn b{color:var(--warn)}
figure{margin:16px 0;background:var(--figbg);border:1px solid var(--line);border-radius:8px;padding:8px;overflow-x:auto}
figure svg.fig{width:100%;height:auto;min-width:640px;display:block}
figcaption{color:#444;font-size:13px;padding:6px 4px 0}
.tablewrap{overflow-x:auto;max-height:560px;overflow-y:auto;border:1px solid var(--line);border-radius:8px;margin:8px 0}
table{border-collapse:collapse;font-size:12.5px;width:100%}
th,td{padding:4px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}
th{position:sticky;top:0;background:var(--card);cursor:pointer;user-select:none}
th.asc::after{content:" \\25B2"}th.desc::after{content:" \\25BC"}
td.na{color:var(--muted);font-style:italic}
input.filter{margin:8px;padding:5px 8px;width:min(320px,90%);border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--fg)}
pre{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px;overflow-x:auto;font-size:12px}
ul.caveats li{margin:6px 0}
"""

HTML_JS = """
document.querySelectorAll('table.sortable').forEach(function(t){
  t.querySelectorAll('th').forEach(function(th,i){
    th.addEventListener('click',function(){
      var asc=!th.classList.contains('asc');
      t.querySelectorAll('th').forEach(function(h){h.classList.remove('asc','desc')});
      th.classList.add(asc?'asc':'desc');
      var rows=Array.from(t.tBodies[0].rows);
      rows.sort(function(a,b){
        var x=a.cells[i],y=b.cells[i];
        var xv=x.dataset.v,yv=y.dataset.v;
        var r=(xv!==undefined&&yv!==undefined)?(parseFloat(xv)-parseFloat(yv)):
              (xv!==undefined?-1:(yv!==undefined?1:x.textContent.localeCompare(y.textContent,undefined,{numeric:true})));
        return asc?r:-r;});
      rows.forEach(function(r){t.tBodies[0].appendChild(r)});
    });
  });
});
document.querySelectorAll('input.filter').forEach(function(inp){
  inp.addEventListener('input',function(){
    var q=inp.value.toLowerCase(),t=document.getElementById(inp.dataset.target);
    Array.from(t.tBodies[0].rows).forEach(function(r){r.style.display=r.textContent.toLowerCase().indexOf(q)>=0?'':'none'});
  });
});
"""


def write_html(outdir, title, ctx):
    fig_by_key = {r["key"]: r for r in ctx["figures"].records}
    fdir = ctx["figures"].dir

    def figure(key):
        r = fig_by_key.get(key)
        if not r:
            return ""
        return (f"<figure>{svg_inline(fdir / (key + '.svg'))}<figcaption><b>{html.escape(r['title'])}.</b> "
                f"{html.escape(r['caption'])}</figcaption></figure>")

    s = ctx["summary"]
    cards = [
        (s["n_samples"], "samples"), (s["n_batches"], "results directories"),
        (s["n_low_depth"], "low-depth samples", s["n_low_depth"] > 0),
    ]
    if s.get("n_tier1") is not None:
        cards += [(s["n_tier1"], "samples with a tier 1 ARG", s["n_tier1"] > 0),
                  (s["n_tier1_genes"], "distinct tier 1 genes")]
    if s.get("n_species") is not None:
        cards.append((s["n_species"], "species detected (any sample)"))
    card_html = "".join(
        f"<div class='card{' warn' if (len(c) > 2 and c[2]) else ''}'><b>{c[0]}</b><span>{html.escape(c[1])}</span></div>"
        for c in cards)

    sections = []
    nav = []

    def section(anchor, heading, body):
        nav.append(f"<a href='#{anchor}'>{html.escape(heading)}</a>")
        sections.append(f"<section id='{anchor}'><h2>{html.escape(heading)}</h2>{body}</section>")

    t = ctx["tables"]
    qc_cols = [c for c in t["samples_qc"][0].columns if c not in ctx["meta_cols"] or c in ctx["key_meta"]]
    section("qc", "Samples and sequencing depth",
            figure("qc_read_fate") + "<h3>Per-sample summary</h3>"
            + df_to_html(t["samples_qc"][0][qc_cols], table_id="t_qc"))
    if "taxonomy_long" in t:
        body = (figure("taxonomy_species_heatmap") + figure("taxonomy_genus_heatmap") + figure("taxonomy_genus_bars")
                + figure("watchlist"))
        if "watchlist" in t:
            wl = t["watchlist"][0]
            body += "<h3>Watchlist (relative abundance %)</h3>" + df_to_html(wl[wl["n_samples_detected"] > 0], table_id="t_wl")
        body += ("<h3>Species prevalence (top 100)</h3>"
                 + df_to_html(t["taxonomy_species_prevalence"][0].head(100), index=True, table_id="t_prev"))
        section("taxonomy", "Taxonomy", body)
    if "arg_long" in t:
        al = t["arg_long"][0]
        keep = ["sample", "tier", "category", "gene", "aro_full_name", "gene_family", "drug_class", "identity",
                "assembly_coverage", "supporting_reads", "rpm", "rpkm", "kmer_blast_discordant"]
        keep += [c for c in ("candidate_hosts", "card_plasmid_associated") if c in al.columns]
        body = (figure("arg_priority_heatmap") + figure("arg_tier1_timeline") + figure("arg_drug_class_heatmap")
                + figure("arg_vs_depth"))
        body += "<h3>Tier 1 (clinical priority) calls</h3>" + df_to_html(al[al["tier"] == 1][keep], table_id="t_t1")
        if "group_arg_detection" in t:
            g = t["group_arg_detection"][0]
            body += (f"<h3>Detection per {html.escape(ctx['group_by'])} (tier 1)</h3>"
                     + df_to_html(g[g["tier"] == 1], table_id="t_grp"))
        t4 = al[al["tier"] == 4][["sample", "aro_full_name", "resistance_snp_status", "resistance_snps_found",
                                  "supporting_reads", "rpkm"]]
        body += ("<h3>Mutation-dependent hits (tier 4)</h3><p class='muted'>Only card_snp_present rows carry a "
                 "mutation CARD lists for resistance; the others are detections of the gene, not of resistance.</p>"
                 + df_to_html(t4, table_id="t_t4"))
        body += "<h3>Tier 2 calls</h3>" + df_to_html(al[al["tier"] == 2][keep], max_rows=1500, table_id="t_t2")
        section("args", "Antimicrobial resistance genes (K-MARVEL)", body)

    caveats = [
        "All counts are read pairs. Reads_report.csv counts single reads (both mates), so its numbers are twice those here.",
        f"Samples with fewer than {ctx['min_pairs']:,} clean read pairs are flagged as low depth. A gene or taxon "
        "not found there is reported as unknown (NA_low_depth), not as absent.",
        "ARG RPKM = supporting read pairs / (reference gene length in kb) / (clean read pairs in millions). "
        "It is comparable between samples, but it is not copies per cell.",
        "Read-level allele calls are unstable among close variants (TEM-1 vs TEM-117), so genes are summarised with "
        "the allele number removed. ESBL status of TEM/SHV depends on the allele and is not called.",
        "Tier 4 (variant-model) hits are detections of a gene that becomes resistant by mutation. Only rows "
        "with resistance_snp_status 'card_snp_present' carry a mutation CARD lists for resistance, and CARD also "
        "lists some common polymorphisms.",
        "Alleles of one gene in one sample usually share their reads, so a gene is represented by its "
        "best-supported allele; allele rows are never summed.",
        "Candidate hosts are species that co-occur in the sample and are known CARD carriers of the ARO. They "
        "are not linkage: confirming the host needs assembly or long reads.",
        "The sylph classified fraction is skewed upward when many reads come from genomes absent from the "
        "database (BugBuster --sylph_read_seq_id corrects it). Taxonomic abundance is relative to the profiled "
        "community.",
        "Tier rules and the watchlist are editable tables (assets/arg_priority_tiers.tsv, assets/summary_watchlist.tsv).",
    ]
    methods = ("<ul class='caveats'>" + "".join(f"<li>{html.escape(c)}</li>" for c in caveats) + "</ul>"
               + "<h3>Inputs</h3><pre>" + html.escape("\n".join(ctx["inputs"])) + "</pre>")
    if ctx["versions"]:
        methods += "<h3>Software versions</h3><pre>" + html.escape(ctx["versions"]) + "</pre>"
    section("methods", "Methods and caveats", methods)

    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title>
<style>{HTML_CSS}</style></head><body>
<header><h1>{html.escape(title)}</h1><div class="muted">BugBuster results summary · generated {date.today().isoformat()}
· {s['n_samples']} samples from {s['n_batches']} results director{'y' if s['n_batches'] == 1 else 'ies'}</div>
<div class="cards">{card_html}</div></header>
<nav>{''.join(nav)}</nav><main>{''.join(sections)}</main><script>{HTML_JS}</script></body></html>"""
    path = outdir / "results_summary.html"
    path.write_text(doc)
    return path


# --------------------------------------------------------------------------- main

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", nargs="+", required=True, help="BugBuster results directories (one per batch)")
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

    batches, seen, labels = [], {}, set()
    for i, r in enumerate(args.results):
        root = Path(r)
        if not root.is_dir():
            die(f"{root} is not a directory")
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
