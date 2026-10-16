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



RGI_ALLELE_COLUMNS = [
    "ARO Term", "ARO Accession", "Reference Model Type", "Reference DB",
    "Resistomes & Variants: Observed Pathogen(s)", "All Mapped Reads", "Percent Coverage",
    "Reference Length", "AMR Gene Family", "Drug Class", "Resistance Mechanism",
]
RGI_KMER_COLUMNS = ["ARO term", "CARD*kmer Prediction"]


def load_rgi(root):
    """RGI bwt allele tables and (optional) RGI k-mer gene tables
    -> (alleles, k-mer predictions, samples rgi bwt ran on, samples rgi kmer_query ran on).

    The allele table lists every reference sequence with at least one mapped read; there is no
    detection call, so the coverage cutoff is applied later (annotate_rgi). A sample where rgi bwt
    found nothing has the low-coverage fallback table (bin/rgi_bwt_lowcov_fallback.sh), whose NA
    row and different header mean "no hits".
    """
    rdir = root / "05_arg_prediction" / "read_level" / "rgi"
    if not rdir.is_dir():
        return None, None, set(), set()
    suffix = "_rgi_bwt.allele_mapping_data.txt"
    frames, ran = [], set()
    for path in sorted(rdir.glob(f"*/*{suffix}")):
        sample = path.name[: -len(suffix)]
        ran.add(sample)
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        values = df.replace("NA", "").astype(str).apply(lambda col: col.str.strip())
        if df.empty or not (values != "").any().any():
            continue
        missing = [c for c in RGI_ALLELE_COLUMNS if c not in df.columns]
        if missing:
            die(f"{path} lacks columns: {', '.join(missing)}")
        df = df[RGI_ALLELE_COLUMNS]
        df = df[df["ARO Accession"].str.fullmatch(r"\d+")]
        df.insert(0, "sample", sample)
        frames.append(df)
    alleles = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["sample"] + RGI_ALLELE_COLUMNS)
    kdir = root / "05_arg_prediction" / "read_level" / "rgi_kmer"
    kframes, kran = [], set()
    for path in sorted(kdir.glob("*/*_rgi_kmer_*mer_analysis.gene.txt")) if kdir.is_dir() else []:
        sample = re.sub(r"_rgi_kmer_\d+mer_analysis\.gene\.txt$", "", path.name)
        kran.add(sample)
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if df.empty:
            continue
        missing = [c for c in RGI_KMER_COLUMNS if c not in df.columns]
        if missing:
            die(f"{path} lacks columns: {', '.join(missing)}")
        df = df[RGI_KMER_COLUMNS].copy()
        df.insert(0, "sample", sample)
        kframes.append(df)
    kmer = pd.concat(kframes, ignore_index=True) if kframes else pd.DataFrame(columns=["sample"] + RGI_KMER_COLUMNS)
    return alleles, kmer, ran, kran


def parse_kmer_prediction(text):
    """'Klebsiella pneumoniae (chromosome or plasmid): 83; Escherichia coli (plasmid): 2; ' -> {taxon: reads}."""
    out = {}
    for part in str(text or "").split(";"):
        m = re.match(r"\s*(.+?)\s*\([^)]*\)\s*:\s*(\d+)\s*$", part)
        if m:
            out[m.group(1)] = out.get(m.group(1), 0) + int(m.group(2))
    return out

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
    arg["gene_rpkm"] = arg["rpkm"].where(arg["gene_representative"])
    arg.insert(1, "method", "kmarvel")
    add_candidate_hosts(arg, prevalence, species_by_sample)
    return arg


def add_candidate_hosts(arg, prevalence, species_by_sample):
    """card_plasmid_associated and candidate_hosts (co-occurring CARD carrier species) columns, in place."""
    if prevalence is None:
        return
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


def watchlist_matcher(watch):
    """Watchlist table -> function(taxon name) -> True when the taxon is a watchlist species or in a watchlist genus."""
    if watch is None or watch.empty:
        return lambda name: False
    species = {norm_taxon(t) for t, r in zip(watch["taxon"], watch["rank"]) if r.lower() == "species"}
    genera = {norm_taxon(t) for t, r in zip(watch["taxon"], watch["rank"]) if r.lower() == "genus"}

    def match(name):
        n = norm_taxon(name)
        return n in species or (n.split(" ")[0] in genera if n else False)
    return match


def host_check(hosts, is_watch, has_kmer):
    """Classify RGI k-mer hosts {taxon: reads} by the top host: watchlist_host, non_watchlist_host,
    genus_level_host (the top prediction names a genus only), no_kmer_host, or no_kmer_data."""
    if not has_kmer:
        return "no_kmer_data"
    if not hosts:
        return "no_kmer_host"
    # the same order rgi_kmer_hosts lists them in: most reads, then name
    top = sorted(hosts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
    if is_watch(top):
        return "watchlist_host"
    return "genus_level_host" if len(top.split(" ")) == 1 else "non_watchlist_host"


def annotate_rgi(alleles, kmer, kmer_ran, samples, rules, prevalence, species_by_sample, watch, min_coverage):
    """RGI bwt alleles -> one row per sample x ARO in the arg_long layout, plus RGI columns.

    coverage (assembly_coverage) = Percent Coverage of the ARO's best-covered reference sequence;
    reads = All Mapped Reads summed over its sequences (KMA assigns each read to one sequence, so,
    unlike K-MARVEL alleles, they do not share reads); read pairs = reads / 2, since RGI counts
    mates separately. detected = coverage >= min_coverage. Genes (alleles collapsed) are detected
    when any of their AROs is; gene_rpkm sums the RPKM of all the gene's AROs, detected or not.
    """
    if alleles is None or alleles.empty:
        return pd.DataFrame()
    a = alleles.copy()
    a["reads"] = pd.to_numeric(a["All Mapped Reads"], errors="coerce").fillna(0)
    a["cov"] = pd.to_numeric(a["Percent Coverage"], errors="coerce").fillna(0)
    a["length"] = pd.to_numeric(a["Reference Length"], errors="coerce")
    a = a.sort_values(["sample", "ARO Accession", "cov", "reads"], ascending=[True, True, False, False])
    best = a.drop_duplicates(["sample", "ARO Accession"]).set_index(["sample", "ARO Accession"])
    reads = a.groupby(["sample", "ARO Accession"])["reads"].sum()
    best = best.assign(reads=reads).reset_index()
    source = best["Reference DB"].where(best["Reference DB"] == "CARD",
                                         "WildCARD: " + best["Resistomes & Variants: Observed Pathogen(s)"])
    arg = pd.DataFrame({
        "sample": best["sample"], "method": "rgi",
        "aro_id": "ARO:" + best["ARO Accession"], "aro_full_name": best["ARO Term"],
        "gene_family": best["AMR Gene Family"], "drug_class": best["Drug Class"],
        "resistance_mechanism": best["Resistance Mechanism"], "model": best["Reference Model Type"],
        "identity": np.nan, "assembly_coverage": best["cov"], "supporting_reads": best["reads"] / 2,
        "median_kmer_depth": np.nan, "gene_length_nt": best["length"],
        "blast_predicted_aro_name": "", "blast_predicted_aro_id": "", "mutation_data": "",
        "rgi_mapped_reads": best["reads"].astype("int64"), "rgi_best_ref": source,
    })
    arg.loc[arg["gene_length_nt"] == 0, "gene_length_nt"] = np.nan
    clean = arg["sample"].map(samples["clean_pairs"]) if "clean_pairs" in samples else np.nan
    arg["rpm"] = arg["supporting_reads"] / clean * 1e6
    arg["rpkm"] = arg["rpm"] / (arg["gene_length_nt"] / 1000)
    tiers = arg.apply(lambda r: assign_tier(r, rules, None), axis=1, result_type="expand")
    arg["tier"], arg["category"], arg["tier_note"] = tiers[0].astype(int), tiers[1], tiers[2]
    arg["gene"] = [gene_name(n, f) for n, f in zip(arg["aro_full_name"], arg["gene_family"])]
    arg["resistance_snp_status"] = np.where(arg["tier"] == 4, "not_checked_rgi", "")
    arg["resistance_snps_found"] = ""
    arg["kmer_blast_discordant"] = ""
    arg["rgi_detected"] = arg["assembly_coverage"] >= min_coverage
    gene_cov = arg.groupby(["sample", "gene"])["assembly_coverage"].transform("max")
    arg["rgi_gene_coverage"] = gene_cov
    arg["rgi_gene_detected"] = gene_cov >= min_coverage
    gene_rpkm = arg.groupby(["sample", "gene"])["rpkm"].transform("sum")
    gene_pairs = arg.groupby(["sample", "gene"])["supporting_reads"].transform("sum")
    # the representative ARO of a gene: detected first, then best covered, then most reads
    arg = arg.sort_values(["sample", "gene", "rgi_detected", "assembly_coverage", "supporting_reads"],
                          ascending=[True, True, False, False, False])
    arg["gene_representative"] = ~arg.duplicated(["sample", "gene"])
    arg["gene_rpkm"] = gene_rpkm.reindex(arg.index).where(arg["gene_representative"])
    arg["rgi_gene_read_pairs"] = gene_pairs.reindex(arg.index)
    # pathogen of origin: RGI k-mer predictions per ARO term, summed over the gene's AROs
    hosts_by = {}
    if not kmer.empty:
        name_gene = dict(zip(zip(arg["sample"], arg["aro_full_name"]), arg["gene"]))
        for s, term, pred in zip(kmer["sample"], kmer["ARO term"], kmer["CARD*kmer Prediction"]):
            g = name_gene.get((s, term))
            if g is None:
                continue
            acc = hosts_by.setdefault((s, g), {})
            for taxon, n in parse_kmer_prediction(pred).items():
                acc[taxon] = acc.get(taxon, 0) + n
    is_watch = watchlist_matcher(watch)
    hosts_txt, checks = [], []
    for s, g in zip(arg["sample"], arg["gene"]):
        h = hosts_by.get((s, g), {})
        top = sorted(h.items(), key=lambda kv: (-kv[1], kv[0]))[:3]
        hosts_txt.append("; ".join(f"{t} ({n})" for t, n in top))
        checks.append(host_check(h, is_watch, s in kmer_ran))
    arg["rgi_kmer_hosts"], arg["rgi_host_check"] = hosts_txt, checks
    add_candidate_hosts(arg, prevalence, species_by_sample)
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
    rpkm = arg.pivot_table(index="gene", columns="sample", values="gene_rpkm", aggfunc="sum")
    rpkm_numeric = rpkm.reindex(columns=cols)
    gene_rpkm = fill_missing(rpkm)
    gene_rpkm.insert(0, "category", gene_meta.loc[gene_rpkm.index, "category"])
    gene_rpkm.insert(0, "tier", gene_meta.loc[gene_rpkm.index, "tier"])
    gene_rpkm = gene_rpkm.sort_values(["tier"], kind="stable")
    out["arg_gene_rpkm"] = gene_rpkm

    fam = arg.pivot_table(index="gene_family", columns="sample", values="gene_rpkm", aggfunc="sum")
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
    dc = counted.pivot_table(index=["tier_group", "drug_class"], columns="sample", values="gene_rpkm", aggfunc="sum")
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
        # typed columns even when no group has a tier 1-2 gene
        out["group_arg_detection"] = pd.DataFrame(rows, columns=[
            group_by, "gene", "tier", "category", "n_detected", "n_evaluable", "n_samples",
            "first_detected", "last_detected", "pattern"])
    return out


METHOD_LABELS = {"kmarvel": "K-MARVEL", "rgi": "RGI"}


def method_combo(methods):
    labels = [METHOD_LABELS.get(m, m) for m in methods]
    return " + ".join(labels) if len(labels) > 1 else f"{labels[0]} only"


def compare_methods(calls, rgi_all, ran_by_method):
    """Read-level ARG methods side by side, per sample x gene (alleles collapsed), on the samples every
    method ran on.

    calls: {method: annotated calls}, each method's detections only; rgi_all: every RGI ARO with reads
    (detected or not), so a gene RGI saw below the coverage cutoff is told apart from one it did not see.
    Returns (per sample x gene table, overlap counts table, {gene: methods} for the UpSet plot).
    """
    methods = [m for m in ("kmarvel", "rgi") if m in calls] + sorted(m for m in calls if m not in ("kmarvel", "rgi"))
    common = set.intersection(*[set(ran_by_method[m]) for m in methods])
    found = {}
    for m in methods:
        d = calls[m]
        if d.empty:
            continue
        d = d[d["gene_representative"] & d["sample"].isin(common)]
        for r in d.to_dict("records"):
            found.setdefault((r["sample"], r["gene"]), {})[m] = r
    rgi_gene = {}
    if rgi_all is not None and not rgi_all.empty:
        rep = rgi_all[rgi_all["gene_representative"] & rgi_all["sample"].isin(common)]
        rgi_gene = {(r["sample"], r["gene"]): r for r in rep.to_dict("records")}
    rows = []
    for (s, g), hits in found.items():
        first = hits[methods[0]] if methods[0] in hits else next(iter(hits.values()))
        row = {"sample": s, "gene": g, "tier": first["tier"], "category": first["category"],
               "gene_family": first["gene_family"], "detected_by": method_combo([m for m in methods if m in hits])}
        for m in methods:
            row[f"detected_{m}"] = m in hits
            pairs = None
            if m in hits:
                pairs = hits[m]["rgi_gene_read_pairs"] if m == "rgi" else hits[m]["supporting_reads"]
            row[f"{m}_read_pairs"] = pairs
        rg = rgi_gene.get((s, g))
        if "rgi" in methods:
            row["rgi_gene_coverage"] = rg["rgi_gene_coverage"] if rg else 0.0
            row["rgi_below_cutoff"] = bool(rg) and "rgi" not in hits
            row["rgi_best_ref"] = rg["rgi_best_ref"] if rg else ""
            row["rgi_kmer_hosts"] = rg["rgi_kmer_hosts"] if rg else ""
            row["rgi_host_check"] = rg["rgi_host_check"] if rg else ""
        rows.append(row)
    comp = pd.DataFrame(rows)
    if comp.empty:
        return comp, pd.DataFrame(), {}
    comp = comp.sort_values(["sample", "tier", "gene"]).reset_index(drop=True)

    genes = {}
    for g, sub in comp.groupby("gene"):
        genes[g] = tuple(m for m in methods if sub[f"detected_{m}"].any())
    gene_meta = comp.drop_duplicates("gene").set_index("gene")[["tier", "category"]]
    below = comp["rgi_below_cutoff"].astype(bool) if "rgi_below_cutoff" in comp else pd.Series(False, comp.index)
    below_det = set(zip(comp.loc[below, "sample"], comp.loc[below, "gene"]))
    below_gene = set(comp.loc[below, "gene"])
    combos = [method_combo(c) for c in _combos(methods)]
    split_below = "rgi" in methods and "kmarvel" in methods
    out = []

    def count(unit, level, members, below_keys):
        groups = {}
        for key, group, combo in members:
            groups.setdefault(group, []).append((key, combo))
        for group in sorted(groups, key=lambda x: (x != "all", x)):
            rec = {"unit": unit, "level": level, "group": group}
            rec.update({c: 0 for c in combos})
            for _, combo in groups[group]:
                rec[combo] += 1
            if split_below:
                only = method_combo(["kmarvel"])
                rec[f"{only}, RGI below cutoff"] = sum(1 for k, c in groups[group] if c == only and k in below_keys)
            rec["total"] = len(groups[group])
            out.append(rec)

    tier = lambda g: f"tier {gene_meta.at[g, 'tier']}"  # noqa: E731
    cat = lambda g: f"T{gene_meta.at[g, 'tier']} {gene_meta.at[g, 'category']}"  # noqa: E731
    for level, keyf in (("all", lambda g: "all"), ("tier", tier), ("category", cat)):
        count("genes", level, [(g, keyf(g), method_combo(ms)) for g, ms in genes.items()], below_gene)
    det = list(zip(zip(comp["sample"], comp["gene"]), comp["detected_by"]))
    for level, keyf in (("all", lambda g: "all"), ("tier", tier)):
        count("detections (sample x gene)", level, [(k, keyf(k[1]), c) for k, c in det], below_det)
    return comp, pd.DataFrame(out), genes


def _combos(methods):
    """Every non-empty subset of the methods, largest first (the column order of the overlap table)."""
    from itertools import combinations
    subsets = []
    for k in range(len(methods), 0, -1):
        subsets += [list(c) for c in combinations(methods, k)]
    return subsets
