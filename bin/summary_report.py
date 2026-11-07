"""The report model of bin/results_summary.py: what the PDF (bin/summary_pdf.py, Typst) and the HTML site
(bin/summary_site.py, Jinja) show, built once so the two differ only in their prose.

One fixed section order; a section is left out when its stage did not run (or found nothing). Each section holds its
facts (numbers the prose cites), its figures and tables, each with a number given here in document order, and its
captions. Tables are formatted to text here: the HTML shows up to `html_rows` rows (sortable, filterable), the PDF
at most PDF_MAX_ROWS rows and PDF_MAX_COLS columns, with a pointer to the full TSV. Methods list every tool and
database of software_versions.yml with the references of assets/report/references.yml.
"""

import csv
from pathlib import Path

import numpy as np
import pandas as pd

from summary_args import METHOD_LABELS
from summary_common import NA_LOW_DEPTH, die, generated_date, warn
from summary_docs import TIER_MEANING, stage_pitfalls, table_doc, tier_categories
from summary_figures import short_count

PDF_MAX_ROWS = 60
PDF_MAX_COLS = 14
PORTRAIT_MM = 183  # text width of a portrait page; a wider PDF table goes on a landscape page

# (id, title, HTML page)
SECTIONS = [
    ("overview", "Overview", "index.html"),
    ("taxonomy", "Taxonomy", "taxonomy.html"),
    ("args", "Antimicrobial resistance genes", "args.html"),
    ("compare", "Comparison of ARG methods", "compare.html"),
    ("sarg", "ARGs-OAP (SARG)", "sarg.html"),
    ("functional", "Functional profile (read level)", "functional.html"),
    ("assembly", "Assembly", "assembly.html"),
    ("mags", "MAGs", "mags.html"),
    ("contig_args", "ARGs on contigs", "contig_args.html"),
    ("bin_args", "ARGs of MAGs", "bin_args.html"),
    ("contig_functional", "Functional profile (contig level)", "contig_functional.html"),
    ("methods", "Methods", "methods.html"),
    ("dictionary", "Data dictionary", "dictionary.html"),
]
PAGE = {sid: page for sid, _, page in SECTIONS}


def fmt(v):
    """A cell as text: floats to 3 significant digits (thousands separated from 1,000), large ints separated,
    missing as blank."""
    if v is None:
        return ""
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (float, np.floating)):
        if np.isnan(v):
            return ""
        if v == 0:
            return "0"
        return f"{v:.3g}" if abs(v) < 1000 else f"{v:,.0f}"
    if isinstance(v, (int, np.integer)):
        return f"{v:,}" if abs(v) >= 10000 else str(v)
    if v is pd.NA or v is pd.NaT:
        return ""
    return str(v)


def table_width_mm(header, rows, pt=6.5, wrap=36):
    """Estimated natural width of a PDF table: per column its longest cell (long text wraps at about wrap
    characters), at pt type, plus the cell padding."""
    width = 0.0
    for i, h in enumerate(header):
        longest = max([len(h)] + [len(r[i]) for r in rows])
        width += min(longest, wrap) * pt * 0.3528 * 0.56 + 2.5
    return width


class Numbering:
    def __init__(self):
        self.fig = 0
        self.tab = 0


def table_entry(nums, df, name, caption, index=False, html_rows=500, pdf=True, pdf_columns=None, tsv=None,
                pdf_rows=PDF_MAX_ROWS):
    """A numbered table: header and text rows for the HTML (up to html_rows) and the PDF (capped)."""
    if df is None:
        return None
    d = df.reset_index() if index else df
    nums.tab += 1
    tsv = tsv if tsv is not None else name
    total = len(d)
    hd = d.head(html_rows)
    out = {"no": nums.tab, "name": name, "caption": caption, "tsv": f"tables/{tsv}.tsv" if tsv else "",
           "header": [str(c) for c in hd.columns], "rows": [[fmt(v) for v in r] for r in hd.itertuples(index=False)],
           "n_rows": total, "html_note": "", "pdf": None}
    if total > html_rows:
        out["html_note"] = (f"First {html_rows:,} of {total:,} rows; the full table is tables/{tsv}.tsv"
                            if tsv else f"First {html_rows:,} of {total:,} rows.")
    if pdf:
        p = d[pdf_columns] if pdf_columns else d
        notes = []
        if p.shape[1] > PDF_MAX_COLS:
            notes.append(f"first {PDF_MAX_COLS} of {p.shape[1]} columns")
            p = p.iloc[:, :PDF_MAX_COLS]
        elif pdf_columns and len(pdf_columns) < d.shape[1]:
            notes.append(f"{len(pdf_columns)} of {d.shape[1]} columns")
        if total > pdf_rows:
            notes.append(f"first {pdf_rows} of {total:,} rows")
            p = p.head(pdf_rows)
        note = ""
        if notes:
            note = ("Shown: " + " and ".join(notes) + (f"; the full table is tables/{tsv}.tsv and the HTML report."
                                                         if tsv else "; the full table is in the HTML report."))
        rows = [[fmt(v) for v in r] for r in p.itertuples(index=False)]
        header = [str(c) for c in p.columns]
        out["pdf"] = {"header": header, "rows": rows, "note": note,
                      "landscape": table_width_mm(header, rows) > PORTRAIT_MM}
    return out


def figure_entry(nums, fw_records, key):
    r = next((x for x in fw_records if x["key"] == key), None)
    if r is None:
        return None
    nums.fig += 1
    return {"no": nums.fig, "key": key, "title": r["title"], "caption": r["caption"], "width_mm": r["width_mm"],
            "height_mm": r["height_mm"], "interactive": r.get("interactive")}


class Section:
    def __init__(self, sid, nums, records):
        self.id = sid
        self.title = dict((s, t) for s, t, _ in SECTIONS)[sid]
        self.page = PAGE[sid]
        self.nums, self.records = nums, records
        self.facts, self.figures, self.tables, self.lists = {}, {}, {}, {}

    def fig(self, key):
        e = figure_entry(self.nums, self.records, key)
        if e:
            self.figures[key] = e
        return e

    def table(self, slot, df, name, caption, **kw):
        if df is None:
            return None
        e = table_entry(self.nums, df, name, caption, **kw)
        self.tables[slot] = e
        return e

    def as_dict(self):
        return {"id": self.id, "title": self.title, "page": self.page, "facts": self.facts, "figures": self.figures,
                "tables": self.tables, "lists": self.lists}


# ---- text shared by both documents (lists, captions, notes): one copy here

def reading_guide(ctx):
    """Units, ARG method notes and tier rows for the methods section."""
    s = ctx["summary"]
    units = ["Counts are read pairs (Reads_report.csv counts single reads, both mates).",
             "K-MARVEL / RGI ARG abundance is RPKM: read pairs per kb of reference gene per million clean read pairs, "
             "comparable between samples of different depth."]
    if ctx.get("sarg"):
        units.append("ARGs-OAP abundance is copies per cell (SARG copies / ARGs-OAP cells), a different unit: do not "
                     "add it to RPKM.")
    if (ctx.get("functional") or {}).get("unit") == "cpge":
        units.append("The read functional profile is in copies per genome equivalent (CPGE: reads per kb / "
                     "MicrobeCensus genome equivalents).")
    units.append(f"Samples with fewer than {ctx['min_pairs']:,} clean read pairs are low depth: a gene or taxon not "
                 f"found there is reported as {NA_LOW_DEPTH} (unknown), not as absent.")
    methods = []
    am = s.get("arg_methods", [])
    if "K-MARVEL" in am:
        methods.append(("K-MARVEL", "assembles resistance genes from the reads and reports a call per gene (CARD), "
                                    "with its supporting read pairs."))
    if "RGI" in am:
        methods.append(("RGI bwt", "maps the reads to CARD and WildCARD and lists every gene with at least one "
                                   "mapped read, fragments of related genes included; this summary counts a gene "
                                   f"when reads cover at least {ctx.get('rgi_min_coverage', 80):g} % of its "
                                   "best-covered reference sequence."))
        methods.append(("RGI k-mer (pathogen of origin)", "assigns a gene's reads to the species they come from "
                        "(rgi_kmer_hosts); rgi_host_check says whether the top species is a watchlist pathogen. "
                        "WildCARD holds resistance genes of gut commensals too, so a tier 1 gene whose reads point "
                        "to a commensal may be a related gene, not the clinical one."))
    if ctx.get("sarg"):
        methods.append(("ARGs-OAP", "searches ARG-like reads against SARG and divides their copies (reads weighted "
                                    "by the share of the reference they cover) by its own cell estimate. SARG names "
                                    "differ from CARD's, so ARGs-OAP has its own section, without tiers."))
    tiers = []
    if s.get("n_tier1") is not None:
        for tier, cats in tier_categories(ctx.get("arg_priority", "")).items():
            label, meaning = TIER_MEANING[tier]
            tiers.append([f"Tier {tier}", label, meaning, "; ".join(cats)])
    return units, methods, tiers


def caveats(ctx):
    out = [
        "All counts are read pairs. Reads_report.csv counts single reads (both mates), so its numbers are twice "
        "those here.",
        f"Samples with fewer than {ctx['min_pairs']:,} clean read pairs are flagged as low depth. A gene or taxon not "
        f"found there is reported as unknown ({NA_LOW_DEPTH}), not as absent.",
        "ARG RPKM = supporting read pairs / (reference gene length in kb) / (clean read pairs in millions). It is "
        "comparable between samples, but it is not copies per cell.",
        "Read-level allele calls are unstable among close variants (TEM-1 vs TEM-117), so genes are summarised with "
        "the allele number removed. ESBL status of TEM/SHV depends on the allele and is not called.",
        "Tier 4 (variant-model) hits are detections of a gene that becomes resistant by mutation. Only rows with "
        "resistance_snp_status 'card_snp_present' carry a mutation CARD lists for resistance, and CARD also lists "
        "some common polymorphisms.",
        "Alleles of one gene in one sample usually share their reads, so a gene is represented by its best-supported "
        "allele; allele rows are never summed.",
        "Candidate hosts are species that co-occur in the sample and are known CARD carriers of the ARO. They are "
        "not linkage: confirming the host needs assembly or long reads.",
        f"RGI bwt maps reads to CARD and WildCARD and reports every gene with at least one mapped read; the summary "
        f"counts an RGI gene when reads cover at least {ctx.get('rgi_min_coverage', 80):g} % of its best-covered "
        "reference sequence (--rgi_summary_min_coverage). RGI counts mates separately, so its read pairs are mapped "
        "reads / 2; a gene's RGI abundance sums its AROs, because KMA assigns each read to one reference.",
        "WildCARD includes resistance-gene sequences found in the genomes of many species, gut commensals among "
        "them, so RGI can report a gene (poxtA, vanD-type, ...) from a related sequence of a commensal. RGI's k-mer "
        "pathogen-of-origin call (rgi_kmer_hosts, from the reads themselves) helps: rgi_host_check says whether the "
        "top host is a watchlist pathogen.",
        "The sylph classified fraction is skewed upward when many reads come from genomes absent from the database "
        "(BugBuster --sylph_read_seq_id corrects it). Taxonomic abundance is relative to the profiled community.",
        "Tier rules and the watchlist are editable tables (assets/arg_priority_tiers.tsv, "
        "assets/summary_watchlist.tsv).",
    ]
    if ctx.get("watchlist") is not None:
        out.append(f"Watchlist: the taxa in the watchlist table, from {ctx.get('watchlist_path', '')}. The default "
                   "list holds the ESKAPE pathogens, the bacteria of the WHO Bacterial Priority Pathogens List 2024 "
                   "and of the CDC Antibiotic Resistance Threats Report 2019, and a few other species of clinical "
                   "interest; each row names its source. A species row matches the species, a genus row every "
                   "species of the genus, GTDB suffixes removed.")
    return out + stage_pitfalls(ctx) + [n for n in ctx.get("method_notes", [])]


FLAG_TITLES = {"tier1": "Tier 1 genes", "single": "Tier 1 genes found by one method only",
               "homolog": "Tier 1 calls that may be commensal relatives", "low_depth": "Depth"}
FLAG_SECTION = {"args": "args", "compare": "compare", "qc": "overview"}


# ---- references

def load_references(assets):
    """references.yml ({key: {author, title, journal, year, volume, pages, doi, type}}) and the tool -> keys map."""
    try:
        import yaml
    except ImportError:
        die("the report's references need pyyaml (conda-forge pyyaml; the pipeline's report image has it)")
    refs_path, map_path = Path(assets) / "references.yml", Path(assets) / "tool_references.tsv"
    refs = yaml.safe_load(refs_path.read_text()) if refs_path.exists() else {}
    tool_refs = {}
    if map_path.exists():
        with open(map_path) as fh:
            for row in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
                tool_refs[row["tool"]] = {"kind": row.get("kind", "tool"),
                                          "keys": [k for k in row.get("references", "").split(";") if k],
                                          "name": row.get("name") or row["tool"]}
    for t, r in tool_refs.items():
        for k in r["keys"]:
            if k not in refs:
                warn(f"tool_references.tsv: {t} cites {k}, which references.yml does not have")
    return refs or {}, tool_refs


def tools_rows(batches, tool_refs):
    """One row per tool / database and version: name, version, kind, processes (and batches when several), keys.
    Helper entries (python, pandas, bash, ...) are left out: tool_references.tsv marks them kind 'helper'."""
    seen = {}
    for b in batches:
        for proc, tools in b["processes"].items():
            for tool, ver in tools.items():
                info = tool_refs.get(tool, {"kind": "tool", "keys": [], "name": tool})
                if info["kind"] == "helper":
                    continue
                key = (tool, ver)
                rec = seen.setdefault(key, {"tool": info["name"], "id": tool, "version": ver, "kind": info["kind"],
                                            "processes": set(), "batches": set(), "keys": info["keys"]})
                rec["processes"].add(proc)
                rec["batches"].add(b["label"])
    rows = sorted(seen.values(), key=lambda r: (r["kind"] != "tool", r["tool"].lower(), str(r["version"])))
    for r in rows:
        r["processes"] = ", ".join(sorted(r["processes"]))
        r["batches"] = ", ".join(sorted(r["batches"]))
    return rows


def format_reference(r):
    """'Authors. Title. Journal volume, pages (year). doi:...' as plain text (the HTML list)."""
    authors = r.get("author", "")
    parts = [a.strip() for a in authors.split(" and ")] if authors else []
    if len(parts) > 5 or (parts and parts[-1] == "others"):
        names = ", ".join(p for p in parts[:3] if p != "others") + " et al."
    else:
        names = ", ".join(parts[:-1]) + (" & " if len(parts) > 1 else "") + (parts[-1] if parts else "")
    venue = r.get("journal") or r.get("publisher") or ""
    vol = f" {r['volume']}" if r.get("volume") else ""
    pages = f", {str(r['pages']).replace('--', '-')}" if r.get("pages") else ""
    doi = f" doi:{r['doi']}" if r.get("doi") else ""
    return f"{names}. {r.get('title', '')}. {venue}{vol}{pages} ({r.get('year', '')}).{doi}"


def bibtex(refs, keys):
    """BibTeX for Typst's bibliography() of the keys cited, in that order."""
    out = []
    for k in keys:
        r = refs[k]
        kind = r.get("type", "article")
        fields = []
        for f in ("author", "title", "journal", "publisher", "year", "volume", "number", "pages", "doi", "url"):
            if r.get(f) not in (None, ""):
                v = str(r[f])
                fields.append(f"  {f} = {{{{{v}}}}}" if f == "title" else f"  {f} = {{{v}}}")
        out.append(f"@{kind}{{{k},\n" + ",\n".join(fields) + "\n}\n")
    return "\n".join(out)


# ---- the model

def build(ctx, assets):
    """The report model (a JSON-ready dict) from results_summary.py's context."""
    t, s = ctx["tables"], ctx["summary"]
    records = ctx["figures"].records
    nums = Numbering()
    secs = []
    present = {r["section"] for r in records}

    # overview
    ov = Section("overview", nums, records)
    qc = t["samples_qc"][0]
    ov.facts = {"n_samples": s["n_samples"], "n_batches": s["n_batches"], "n_low": s["n_low_depth"],
                "min_pairs": f"{ctx['min_pairs']:,}", "min_pairs_short": short_count(ctx["min_pairs"]),
                "low_samples": ", ".join(s.get("low_depth_samples", [])),
                "batches": ", ".join(b["label"] for b in ctx["batches"])}
    if "clean_pairs" in qc and qc["clean_pairs"].notna().any():
        ov.facts["median_clean"] = f"{int(qc['clean_pairs'].median()):,}"
    cards = [(s["n_samples"], "samples", False),
             (s["n_batches"], "results directories", False),
             (s["n_low_depth"], f"low-depth samples (< {short_count(ctx['min_pairs'])} clean pairs)",
              s["n_low_depth"] > 0)]
    for mark, n in s.get("depth_marks", []):
        cards.append((n, f"samples with >= {short_count(mark)} clean pairs", False))
    if s.get("n_tier1") is not None:
        cards += [(s["n_tier1"], f"samples with a tier 1 ARG ({s.get('arg_primary', '')})", s["n_tier1"] > 0),
                  (s["n_tier1_genes"], f"distinct tier 1 genes ({s.get('arg_primary', '')})", False)]
    if s.get("n_tier1_genes_all_methods") is not None:
        cards.append((s["n_tier1_genes_all_methods"], "tier 1 genes found by every ARG method", False))
    if s.get("n_species") is not None:
        cards.append((s["n_species"], "species detected (any sample)", False))
    ov.lists["cards"] = [{"value": fmt(v), "label": lab, "warn": bool(w)} for v, lab, w in cards]
    att = []
    for kind in ("tier1", "single", "homolog", "low_depth"):
        sub = [f for f in s.get("attention", []) if f[0] == kind]
        if sub:
            att.append({"kind": kind, "title": FLAG_TITLES[kind],
                        "items": [{"text": text, "page": PAGE[FLAG_SECTION.get(a, "overview")]} for _, text, a in sub]})
    ov.lists["attention"] = att
    labels = list(ctx["stage_status"])
    stage_rows = []
    for i, (st, _) in enumerate(ctx["stage_status"][labels[0]]):
        stage_rows.append([st["title"]] + [ctx["stage_status"][b][i][1] for b in labels] + [st["summary"]])
    ov.facts["stages_ran"] = sum(1 for r in stage_rows if any(x.startswith("ran") for x in r[1:-1]))
    ov.facts["stages_total"] = len(stage_rows)
    ov.table("stages", pd.DataFrame(stage_rows, columns=["Stage"] + labels + ["In this summary"]), "stages",
             "Pipeline stages: whether each ran in each results directory (from software_versions.yml) and what "
             "this summary takes from it.", tsv="")
    ov.fig("qc_read_fate")
    qc_cols = [c for c in qc.columns if c not in ctx["meta_cols"] or c in ctx["key_meta"]]
    pdf_qc = [c for c in ["sample", "batch"] + ctx["key_meta"] + ["raw_pairs", "clean_pairs", "pct_clean",
                                                                  "low_depth", "depth_class", "classified_fraction",
                                                                  "n_species", "shannon", "acquired_genes"]
              if c in qc_cols]
    ov.table("samples", qc[qc_cols], "samples_qc", "Per-sample summary: read pairs per step, depth class, profile "
             "and ARG counts (samples_qc; every column is described in the data dictionary).", pdf_columns=pdf_qc)
    secs.append(ov)

    if "taxonomy_long" in t:
        tx = Section("taxonomy", nums, records)
        tx.facts = {"db": str(qc["taxonomy_db"].iloc[0]) if "taxonomy_db" in qc else "",
                    "n_species": s.get("n_species")}
        for k in ("taxonomy_species_heatmap", "taxonomy_genus_heatmap", "taxonomy_genus_bars", "watchlist"):
            tx.fig(k)
        if "watchlist" in t:
            wl = t["watchlist"][0]
            tx.table("watchlist", wl[wl["n_samples_detected"] > 0], "watchlist",
                     "Watchlist taxa detected in at least one sample (relative abundance, %); the full list and "
                     "each taxon's source are in Methods.")
        tx.table("prevalence", t["taxonomy_species_prevalence"][0].head(100), "taxonomy_species_prevalence",
                 "Species prevalence, the 100 most prevalent species.", index=True)
        secs.append(tx)

    if "arg_long" in t:
        ar = Section("args", nums, records)
        primary = s.get("arg_primary", "")
        others = [m for m in s.get("arg_methods", []) if m != primary]
        ar.title = f"Antimicrobial resistance genes ({primary})"
        ar.facts = {"primary": primary, "others": ", ".join(others), "rgi_cov": f"{ctx['rgi_min_coverage']:g}",
                    "n_tier1_samples": s.get("n_tier1"), "n_tier1_genes": s.get("n_tier1_genes"),
                    "comparison": ctx.get("comparison") is not None}
        for k in ("arg_priority_heatmap", "arg_drug_class_heatmap", "arg_vs_depth"):
            ar.fig(k)
        al = t["arg_long"][0]
        al = al[al["method"].map(lambda m: METHOD_LABELS.get(m, m)) == primary]
        t1 = ctx.get("tier1_view")
        if t1 is not None and not t1.empty:
            ar.table("tier1", t1, "arg_long", "Tier 1 (clinical priority) calls with the other method's evidence. "
                     "rgi_host_check: watchlist_host = RGI's k-mer pathogen-of-origin call is a watchlist pathogen; "
                     "non_watchlist_host = another species (often a related gene of a commensal); genus_level_host = "
                     "the reads only resolve a genus.", tsv="arg_long")
        if "group_arg_detection" in t and not t["group_arg_detection"][0].empty:
            g = t["group_arg_detection"][0]
            ar.table("group", g[g["tier"] == 1], "group_arg_detection",
                     f"Tier 1 detection per {ctx['group_by']}.")
        t4 = al[al["tier"] == 4][["sample", "aro_full_name", "resistance_snp_status", "resistance_snps_found",
                                  "supporting_reads", "rpkm"]]
        ar.table("tier4", t4, "arg_long", "Mutation-dependent hits (tier 4). Only card_snp_present rows carry a "
                 "mutation CARD lists for resistance; the others are detections of the gene, not of resistance.")
        keep = ["sample", "tier", "category", "gene", "aro_full_name", "gene_family", "drug_class", "identity",
                "assembly_coverage", "supporting_reads", "rpm", "rpkm", "kmer_blast_discordant"]
        keep += [c for c in ("candidate_hosts", "card_plasmid_associated", "rgi_kmer_hosts", "rgi_host_check")
                 if c in al.columns and al[c].astype(str).str.len().gt(0).any()]
        ar.table("tier2", al[al["tier"] == 2][keep], "arg_long", f"Tier 2 (other acquired) calls of {primary}.",
                 html_rows=1500, pdf_columns=["sample", "gene", "category", "supporting_reads", "rpkm"])
        secs.append(ar)

    comp = ctx.get("comparison")
    if comp is not None:
        co = Section("compare", nums, records)
        co.facts = {"rgi_cov": f"{ctx['rgi_min_coverage']:g}",
                    "methods": " + ".join(METHOD_LABELS.get(m, m) for m in comp["methods"])}
        co.fig("arg_method_upset")
        co.fig("arg_method_by_category")
        co.table("overlap", comp["overlap"], "arg_method_overlap", "Distinct genes by the methods that found them.")
        co.table("comparison", comp["comp"], "arg_method_comparison", "Every gene per sample and method.",
                 html_rows=2000, pdf_columns=[c for c in ("sample", "gene", "tier", "category", "detected_by",
                                                          "kmarvel_read_pairs", "rgi_read_pairs", "rgi_gene_coverage")
                                              if c in comp["comp"].columns])
        secs.append(co)

    if ctx.get("sarg"):
        sg = Section("sarg", nums, records)
        g = ctx["sarg"]
        sg.facts = {"n_samples": g["n_samples"], "n_types": g["n_types"], "n_subtypes": g["n_subtypes"],
                    "not_normalised": ", ".join(g["not_normalised"])}
        sg.fig("arg_sarg_type_heatmap")
        if "arg_sarg_type_cpc" in t:
            sg.table("types", t["arg_sarg_type_cpc"][0], "arg_sarg_type_cpc", "Copies per cell by SARG type.",
                     index=True)
        secs.append(sg)

    if ctx.get("functional"):
        fn = Section("functional", nums, records)
        f = ctx["functional"]
        fn.facts = {"backend": f["backend"], "n_samples": f["n_samples"], "ontologies": ", ".join(f["ontologies"]),
                    "unit": "copies per genome equivalent (CPGE)" if f["unit"] == "cpge" else f["unit"]}
        fn.fig("function_cog_heatmap")
        cols = [c for c in ("sample", "func_reads_assigned", "func_fraction_annotated", "mc_avg_genome_size_bp",
                            "mc_genome_equivalents", "args_oap_ncell") if c in qc.columns]
        fn.table("per_sample", qc[cols], "samples_qc", "Read functional profile per sample, with MicrobeCensus.")
        if "function_read_ko" in t:
            fn.table("ko", t["function_read_ko"][0].head(50), "function_read_ko",
                     "The 50 most abundant KEGG orthologs (read level).", index=True)
        secs.append(fn)

    if "assembly_stats" in t:
        asm = Section("assembly", nums, records)
        asm.facts = {"n_assemblies": len(t["assembly_stats"][0])}
        asm.table("stats", t["assembly_stats"][0], "assembly_stats", "MEGAHIT assemblies: statistics of all "
                  "contigs, and of those kept by the length filter that binning, contig taxonomy, ARGs and functions "
                  "use.")
        if "assembly_read_use" in t:
            asm.table("read_use", t["assembly_read_use"][0], "assembly_read_use", "Reads on the assembly: the share "
                      "of each sample's reads on the contigs and on predicted genes; reads not on contigs are "
                      "outside every contig-based table.")
        secs.append(asm)

    if "mags" in t:
        mg = Section("mags", nums, records)
        m = t["mags"][0]
        tiers = m["quality_tier"].value_counts()
        mg.facts = {"n_mags": len(m), "high": int(tiers.get("high", 0)), "medium": int(tiers.get("medium", 0)),
                    "low": int(tiers.get("low", 0))}
        mg.fig("mag_quality")
        mg.fig("mag_abundance_heatmap")
        mg.table("mags", m, "mags", "Refined MAGs: CheckM2 quality, MIMAG tier (tier_basis: whether Bakta's rRNA and "
                 "tRNA genes were checked), GTDB-Tk lineage and Bakta counts.",
                 pdf_columns=[c for c in ("assembly", "bin", "completeness_pct", "contamination_pct", "quality_tier",
                                          "genome_size_mb", "contigs", "phylum", "genus", "species") if c in m.columns])
        if "mag_abundance" in t:
            mg.table("abundance", t["mag_abundance"][0], "mag_abundance", "MAG depth per sample; pct_sample_bases "
                     "compares samples of one assembly mode, raw depth does not.", html_rows=1000)
        secs.append(mg)

    if "arg_contig_long" in t or "arg_contig_mobile" in t:
        ca = Section("contig_args", nums, records)
        cl = t["arg_contig_long"][0] if "arg_contig_long" in t else None
        ca.facts = {"n_genes": len(cl) if cl is not None else 0,
                    "n_assemblies": int(cl["assembly"].nunique()) if cl is not None else 0}
        ca.fig("arg_contig_hosts")
        if "arg_contig_hosts" in t:
            ca.table("hosts", t["arg_contig_hosts"][0], "arg_contig_hosts", "ARG genes per host taxon, ARG class and "
                     "assembly.")
        if "arg_contig_mobile" in t:
            ca.table("mobile", t["arg_contig_mobile"][0], "arg_contig_mobile", "ARG genes on plasmids, viruses and "
                     "proviruses (geNomad) per assembly and ARG class.")
        if cl is not None:
            ca.table("long", cl, "arg_contig_long", "Every ARG gene on a contig, with the contig's taxon, bin and "
                     "geNomad class.", html_rows=1000,
                     pdf_columns=[c for c in ("assembly", "contig", "arg_gene", "arg_class", "identity_pct",
                                              "host_rank", "host_taxon", "bin", "mobile_class") if c in cl.columns])
        secs.append(ca)

    if "arg_mag_long" in t:
        ba = Section("bin_args", nums, records)
        ba.facts = {"n_genes": len(t["arg_mag_long"][0]), "n_shared": len(t["arg_cluster_sharing"][0])}
        ba.table("sharing", t["arg_cluster_sharing"][0], "arg_cluster_sharing", "ARG clusters (90 % identity) found "
                 "in more than one MAG.")
        ba.table("long", t["arg_mag_long"][0], "arg_mag_long", "ARG genes of the MAGs and their clusters at 90-100 % "
                 "identity.", html_rows=1000)
        secs.append(ba)

    if ctx.get("contig_functional"):
        cf = Section("contig_functional", nums, records)
        c = ctx["contig_functional"]
        cf.facts = {"n_samples": c["n_samples"], "ontologies": ", ".join(c["ontologies"]),
                    "unit": "copies per genome equivalent (CPGE)" if c["unit"] == "cpge" else "TPM"}
        cf.fig("function_contig_cog_heatmap")
        cols = [x for x in ("sample", "pct_reads_on_genes", "func_contig_genes_annotated_pct",
                            "func_contig_annotated_by_abundance_pct") if x in qc.columns]
        cf.table("per_sample", qc[cols], "samples_qc", "Contig functional profile per sample.")
        for ont, title in (("cog", "COG categories (contig level)."),
                           ("ko", "The 50 most abundant KEGG orthologs (contig level).")):
            if f"function_contig_{ont}" in t:
                cf.table(ont, t[f"function_contig_{ont}"][0].head(50), f"function_contig_{ont}", title, index=True)
        secs.append(cf)

    # methods
    me = Section("methods", nums, records)
    units, mnotes, tier_rows = reading_guide(ctx)
    me.lists["units"] = units
    me.lists["arg_methods"] = [{"term": a, "text": b} for a, b in mnotes]
    me.lists["caveats"] = caveats(ctx)
    me.lists["inputs"] = ctx["inputs"]
    me.facts = {"arg_priority": ctx.get("arg_priority_name", ""), "primary": s.get("arg_primary", ""),
                "n_methods": len(s.get("arg_methods", []))}
    if tier_rows:
        me.table("tiers", pd.DataFrame(tier_rows, columns=["Tier", "Name", "Meaning", "Categories (rule notes)"]),
                 "tiers", f"ARG reporting tiers, from {ctx.get('arg_priority_name', '')}: rules match CARD gene names "
                 "or families, top to bottom, first match wins.", tsv="")
    refs, tool_refs = load_references(assets)
    rows = tools_rows(ctx["batches"], tool_refs)
    cited = []
    for r in rows:
        for k in r["keys"]:
            if k in refs and k not in cited:
                cited.append(k)
    multi = len(ctx["batches"]) > 1
    tool_df = pd.DataFrame([{"Tool": r["tool"], "Version": r["version"], "Kind": r["kind"],
                             "Processes": r["processes"], **({"Batches": r["batches"]} if multi else {}),
                             "Ref.": ",".join(str(cited.index(k) + 1) for k in r["keys"] if k in cited)}
                            for r in rows])
    if not tool_df.empty:
        e = me.table("tools", tool_df, "tools", "Software and databases with their versions and the processes that "
                     "used them, from software_versions.yml.", tsv="", html_rows=100000,
                     pdf_rows=100000)
        e["keys"] = [[k for k in r["keys"] if k in cited] for r in rows]
        if e.get("pdf"):
            e["pdf"]["landscape"] = False
    me.lists["references"] = [{"no": i + 1, "key": k, "text": format_reference(refs[k])} for i, k in enumerate(cited)]
    me.facts["bibtex"] = bibtex(refs, cited)
    if ctx.get("watchlist") is not None:
        me.table("watchlist", ctx["watchlist"], "watchlist_definition", f"Watchlist taxa ({ctx.get('watchlist_path')})"
                 ".", tsv="")
    secs.append(me)

    # data dictionary: every table written, with its columns
    dd = Section("dictionary", nums, records)
    entries = []
    for name in t:
        doc = table_doc(name, ctx.get("group_by"))
        if doc is None:
            continue
        unit = ctx.get("table_units", {}).get(name)
        entries.append({"name": name, "title": doc["title"] + (f" ({unit})" if unit else ""),
                        "description": doc["description"],
                        "columns": [{"name": c, "text": d} for c, d in doc["columns"]]})
    dd.lists["tables"] = entries
    secs.append(dd)

    unused = present - {sec.id for sec in secs if sec.figures}
    if unused:
        warn(f"figures of sections not in the report: {sorted(unused)}")
    meta = {"title": ctx["title"], "generated": generated_date(), "n_samples": s["n_samples"],
            "batches": [b["label"] for b in ctx["batches"]],
            "pipeline_version": pipeline_version(ctx["batches"])}
    return {"run": meta, "sections": [sec.as_dict() for sec in secs]}


def pipeline_version(batches):
    vers = sorted({str(v) for b in batches for p, tools in b["processes"].items() for k, v in tools.items()
                   if k in ("pipeline", "bugbuster", "BugBuster")})
    return ", ".join(vers)
