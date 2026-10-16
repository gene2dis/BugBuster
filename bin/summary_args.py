"""ARG loading, tier assignment and ARG matrices for bin/results_summary.py."""

import re
from pathlib import Path

import numpy as np
import pandas as pd

from summary_common import NA_LOW_DEPTH, NA_NOT_RUN, die, norm_taxon, warn


KMARVEL_COLUMNS = [
    "aro_id", "aro_full_name", "gene_family", "drug_class", "resistance_mechanism",
    "model", "identity", "assembly_coverage", "supporting_reads", "median_kmer_depth",
    "blast_predicted_aro_name", "blast_predicted_aro_id", "mutation_data",
    "reference_sequence",
]


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
