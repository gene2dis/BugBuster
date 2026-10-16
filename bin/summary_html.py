"""Table, Excel and HTML output of bin/results_summary.py."""

import html
import re
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from summary_common import warn
from summary_figures import short_count


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
        ("acquired_genes, acquired_rpkm", "distinct tier 1-2 genes and their summed RPKM (primary ARG method)"),
        ("depth_class", "clean read pairs binned at --min-clean-pairs and --depth-marks (e.g. <1M, 1M-10M, >=20M)"),
        ("n_genes_kmarvel, n_genes_rgi, n_genes_all_methods", "distinct genes per read-level ARG method and found "
                                                              "by all of them (when more than one ran)"),
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
        ("method", "kmarvel or rgi; every method's calls are listed, the matrices use the primary method only"),
        ("assembly_coverage", "K-MARVEL: coverage of the assembled gene; RGI: % of the ARO's best-covered reference "
                              "sequence covered by reads"),
        ("gene_rpkm", "RPKM of the gene in the sample, on its gene_representative row (K-MARVEL: the best-supported "
                      "allele; RGI: summed over the gene's AROs)"),
        ("rgi_mapped_reads, rgi_best_ref", "RGI: mapped reads (mates counted separately) and the source of the best "
                                           "reference sequence (CARD, or WildCARD and the species it came from)"),
        ("rgi_gene_coverage", "RGI: best coverage of any ARO of the gene"),
        ("rgi_kmer_hosts, rgi_host_check", "RGI k-mer pathogen of origin (top hosts, reads) and whether the top host "
                                           "is a watchlist taxon (watchlist_host, non_watchlist_host, "
                                           "genus_level_host, no_kmer_host, no_kmer_data)"),
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
        ("arg_method_comparison", "per sample and gene (alleles collapsed): detected by each method, read pairs, RGI "
                                  "coverage, rgi_below_cutoff (RGI reads with partial coverage), RGI hosts"),
        ("arg_method_overlap", "distinct genes, and sample x gene detections, found by every method or by one only, "
                               "overall, per tier and per category"),
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
ul.tight{margin:0;padding-left:16px}ul.tight li{margin:1px 0}
#guide td{white-space:normal;vertical-align:top}
details.attention{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:8px 14px;margin:4px 0 8px}
details.attention summary{font-weight:600;cursor:pointer}
details.attention>ul{margin:6px 0;padding-left:18px}details.attention ul ul{padding-left:18px;font-size:13.5px}
details.attention li.homolog>b,details.attention li.single>b{color:var(--warn)}
details.attention a{color:inherit;text-decoration:none}details.attention a:hover{text-decoration:underline}
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


TIER_MEANING = {
    1: ("Clinical priority", "last-line or highly transmissible resistance, usually on mobile elements"),
    2: ("Other acquired", "acquired genes of common resistance (ESBL, AmpC, aminoglycoside-modifying enzymes, "
                          "tetracycline, macrolide, sulfonamide/trimethoprim genes, ...)"),
    3: ("Intrinsic, efflux or regulator", "chromosomal or housekeeping genes that most strains of a species carry; "
                                          "they say little about acquired resistance"),
    4: ("Mutation dependent", "genes that confer resistance only with a specific mutation (gyrA, parC, rpoB, ...)"),
}


def tier_categories(rules_path):
    """{tier: ['category (note)', ...]} from the tier rule file, in rule order, plus the built-in rules."""
    out = {1: [], 2: [], 3: [], 4: ["mutation_dependent (CARD variant models and AROs with CARD SNP definitions)"]}
    try:
        rules = pd.read_csv(rules_path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    except (OSError, ValueError):
        rules = pd.DataFrame(columns=["tier", "category", "note"])
    seen = set()
    for tier, cat, note in zip(rules.get("tier", []), rules.get("category", []), rules.get("note", [])):
        key = (int(tier), cat)
        if key in seen:
            continue
        seen.add(key)
        out.setdefault(int(tier), []).append(f"{cat} ({note})" if note else cat)
    out[3].append("efflux_permeability (efflux pumps, porins and overexpression models without a rule)")
    out[2].append("acquired_unreviewed (no rule matched)")
    return out


def reading_guide(ctx):
    """The 'Reading this report' section: units, ARG tiers and the read-level ARG methods."""
    s = ctx["summary"]
    rows = []
    for tier, cats in tier_categories(ctx.get("arg_priority", "")).items():
        label, meaning = TIER_MEANING[tier]
        rows.append(f"<tr><td><b>Tier {tier}</b><br>{html.escape(label)}</td><td>{html.escape(meaning)}</td>"
                    f"<td><ul class='tight'>{''.join(f'<li>{html.escape(c)}</li>' for c in cats)}</ul></td></tr>")
    tiers = ("<h3>ARG reporting tiers</h3><p>Every ARG call gets a tier from an editable rule table ("
             + html.escape(Path(ctx.get("arg_priority", "")).name or "arg_priority_tiers.tsv")
             + "): rules match CARD gene names or gene families, top to bottom, first match wins. Tier 1 is what to "
             "look at first; tier 2 is the common acquired resistome; tiers 3 and 4 are mostly background.</p>"
             "<div class='tablewrap'><table><thead><tr><th>Tier</th><th>Meaning</th><th>Categories (rule notes)</th>"
             "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>")
    methods = s.get("arg_methods", [])
    m = []
    if "K-MARVEL" in methods:
        m.append("<li><b>K-MARVEL</b> assembles resistance genes from the reads and reports a call per gene (CARD), "
                 "with its supporting read pairs.</li>")
    if "RGI" in methods:
        m.append(f"<li><b>RGI bwt</b> maps the reads to CARD and WildCARD and lists every gene with at least one "
                 f"mapped read, fragments of related genes included. This summary counts an RGI gene when reads cover "
                 f"at least {ctx.get('rgi_min_coverage', 80):g}% of its best-covered reference sequence.</li>"
                 "<li><b>RGI k-mer (pathogen of origin)</b> assigns a gene's reads to the species they come from "
                 "(<i>rgi_kmer_hosts</i>). <i>rgi_host_check</i> says whether the top species is a watchlist pathogen; "
                 "WildCARD holds resistance genes of gut commensals too, so a tier 1 gene whose reads point to a "
                 "commensal may be a related gene, not the clinical one.</li>")
    method_html = ""
    if m:
        primary = s.get("arg_primary", "")
        method_html = ("<h3>ARG methods</h3><ul>" + "".join(m) + "</ul><p>The ARG figures and matrices use "
                       f"<b>{html.escape(primary)}</b>"
                       + ("; the other method's calls are in arg_long.tsv and in the method comparison section."
                          if len(methods) > 1 else ".") + "</p>")
    units = ("<h3>Units</h3><ul><li>Counts are read pairs (Reads_report.csv counts single reads).</li>"
             "<li>ARG abundance is RPKM: read pairs per kb of reference gene per million clean read pairs, comparable "
             "between samples of different depth.</li>"
             f"<li>Samples with fewer than {ctx['min_pairs']:,} clean read pairs are low depth: a gene or taxon not "
             "found there is reported as NA_low_depth (unknown), not as absent.</li></ul>")
    return units + (tiers if s.get("n_tier1") is not None else "") + method_html


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
        (s["n_low_depth"], f"low-depth samples (<{short_count(ctx['min_pairs'])} clean pairs)", s["n_low_depth"] > 0),
    ]
    for mark, n in s.get("depth_marks", []):
        cards.append((n, f"samples with >={short_count(mark)} clean pairs"))
    if s.get("n_tier1") is not None:
        cards += [(s["n_tier1"], f"samples with a tier 1 ARG ({s.get('arg_primary', '')})", s["n_tier1"] > 0),
                  (s["n_tier1_genes"], "distinct tier 1 genes")]
    if s.get("n_tier1_genes_all_methods") is not None:
        cards.append((s["n_tier1_genes_all_methods"], "tier 1 genes found by every ARG method"))
    if s.get("n_species") is not None:
        cards.append((s["n_species"], "species detected (any sample)"))
    card_html = "".join(
        f"<div class='card{' warn' if (len(c) > 2 and c[2]) else ''}'><b>{c[0]}</b><span>{html.escape(c[1])}</span></div>"
        for c in cards)
    flag_titles = {"tier1": "Tier 1 genes", "single": "Tier 1 genes found by one method only",
                   "homolog": "Tier 1 calls that may be commensal relatives", "low_depth": "Depth"}
    attention = ""
    if s.get("attention"):
        items = []
        for kind in ("tier1", "single", "homolog", "low_depth"):
            sub = [f for f in s["attention"] if f[0] == kind]
            if sub:
                items.append(f"<li class='{kind}'><b>{html.escape(flag_titles[kind])}</b><ul>"
                             + "".join(f"<li><a href='#{a}'>{html.escape(text)}</a></li>" for _, text, a in sub)
                             + "</ul></li>")
        attention = ("<details class='attention' open><summary>Attention</summary><p class='muted'>Tiers, ARG "
                     "methods and the RGI host are explained in <a href='#guide'>Reading this report</a>.</p><ul>"
                     + "".join(items) + "</ul></details>")

    sections = []
    nav = []

    def section(anchor, heading, body):
        nav.append(f"<a href='#{anchor}'>{html.escape(heading)}</a>")
        sections.append(f"<section id='{anchor}'><h2>{html.escape(heading)}</h2>{body}</section>")

    t = ctx["tables"]
    section("guide", "Reading this report", reading_guide(ctx))
    qc_cols = [c for c in t["samples_qc"][0].columns if c not in ctx["meta_cols"] or c in ctx["key_meta"]]
    section("qc", "Samples and sequencing depth",
            figure("qc_read_fate") + "<h3>Per-sample summary</h3>"
            + df_to_html(t["samples_qc"][0][qc_cols], table_id="t_qc"))
    if "taxonomy_long" in t:
        body = (figure("taxonomy_species_heatmap") + figure("taxonomy_genus_heatmap") + figure("taxonomy_genus_bars")
                + figure("watchlist"))
        if "watchlist" in t:
            wl = t["watchlist"][0]
            body += ("<h3>Watchlist (relative abundance %)</h3><p class='muted'>Watchlist taxa detected in at least "
                     "one sample; the full list and the source of each taxon are under <a href='#methods'>Methods</a>.</p>"
                     + df_to_html(wl[wl["n_samples_detected"] > 0], table_id="t_wl"))
        body += ("<h3>Species prevalence (top 100)</h3>"
                 + df_to_html(t["taxonomy_species_prevalence"][0].head(100), index=True, table_id="t_prev"))
        section("taxonomy", "Taxonomy", body)
    if "arg_long" in t:
        al = t["arg_long"][0]
        primary = s.get("arg_primary", "")
        al = al[al["method"].map(lambda m: {"kmarvel": "K-MARVEL", "rgi": "RGI"}.get(m, m)) == primary]
        keep = ["sample", "tier", "category", "gene", "aro_full_name", "gene_family", "drug_class", "identity",
                "assembly_coverage", "supporting_reads", "rpm", "rpkm", "kmer_blast_discordant"]
        keep += [c for c in ("candidate_hosts", "card_plasmid_associated", "rgi_kmer_hosts", "rgi_host_check")
                 if c in al.columns and al[c].astype(str).str.len().gt(0).any()]
        intro = (f"<p class='muted'>Matrices, figures and counts use {html.escape(primary)}"
                 + (f"; arg_long.tsv also holds the calls of {html.escape(', '.join(m for m in s.get('arg_methods', []) if m != primary))}"
                    if len(s.get("arg_methods", [])) > 1 else "")
                 + f". RGI genes are counted at >= {ctx['rgi_min_coverage']:g}% coverage of their best reference "
                 "sequence (RGI itself reports every gene with a mapped read).</p>")
        body = (intro + figure("arg_priority_heatmap") + figure("arg_drug_class_heatmap")
                + figure("arg_vs_depth"))
        if ctx.get("tier1_view") is not None:
            body += ("<h3>Tier 1 (clinical priority) calls</h3><p class='muted'>rgi_host_check: watchlist_host = RGI's "
                     "k-mer pathogen-of-origin call is a watchlist pathogen; non_watchlist_host = another species "
                     "(often a related gene of a commensal); genus_level_host = the reads only resolve a genus.</p>"
                     + df_to_html(ctx["tier1_view"], table_id="t_t1"))
        if "group_arg_detection" in t and not t["group_arg_detection"][0].empty:
            g = t["group_arg_detection"][0]
            body += (f"<h3>Detection per {html.escape(ctx['group_by'])} (tier 1)</h3>"
                     + df_to_html(g[g["tier"] == 1], table_id="t_grp"))
        t4 = al[al["tier"] == 4][["sample", "aro_full_name", "resistance_snp_status", "resistance_snps_found",
                                  "supporting_reads", "rpkm"]]
        body += ("<h3>Mutation-dependent hits (tier 4)</h3><p class='muted'>Only card_snp_present rows carry a "
                 "mutation CARD lists for resistance; the others are detections of the gene, not of resistance.</p>"
                 + df_to_html(t4, table_id="t_t4"))
        body += "<h3>Tier 2 calls</h3>" + df_to_html(al[al["tier"] == 2][keep], max_rows=1500, table_id="t_t2")
        section("args", f"Antimicrobial resistance genes ({primary})", body)
    comp = ctx.get("comparison")
    if comp is not None:
        ov = comp["overlap"]
        body = ("<p class='muted'>Read-level methods on the same samples, gene by gene (allele numbers removed, so "
                "KPC-2 and KPC-97 are both KPC). K-MARVEL reports calls; RGI bwt reports every gene with a mapped "
                f"read, so an RGI gene counts when reads cover >= {ctx['rgi_min_coverage']:g}% of its best reference "
                "sequence. 'RGI below cutoff' marks K-MARVEL genes RGI saw with partial coverage, which is not absence."
                " Tier 4 (mutation-dependent) genes are mostly K-MARVEL calls: RGI bwt reports homolog matches and does "
                "not test resistance mutations. KARGA, KARGVA and ARGs-OAP use other databases (MEGARes, SARG) and are "
                "not compared.</p>"
                + figure("arg_method_upset") + figure("arg_method_by_category")
                + "<h3>Overlap counts</h3>" + df_to_html(ov, table_id="t_ov")
                + "<h3>Per sample and gene</h3>" + df_to_html(comp["comp"], max_rows=2000, table_id="t_cmp"))
        section("compare", "Comparison of ARG methods", body)

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
        f"RGI bwt maps reads to CARD and WildCARD and reports every gene with at least one mapped read; the summary "
        f"counts an RGI gene when reads cover at least {ctx.get('rgi_min_coverage', 80):g}% of its best-covered "
        "reference sequence (--rgi_summary_min_coverage). RGI counts mates separately, so its read pairs are mapped "
        "reads / 2; a gene's RGI abundance sums its AROs, because KMA assigns each read to one reference.",
        "WildCARD includes resistance-gene sequences found in the genomes of many species, gut commensals among them, "
        "so RGI can report a gene (poxtA, vanD-type, ...) from a related sequence of a commensal. RGI's k-mer "
        "pathogen-of-origin call (rgi_kmer_hosts, from the reads themselves) helps: rgi_host_check says whether the "
        "top host is a watchlist pathogen.",
        "The sylph classified fraction is skewed upward when many reads come from genomes absent from the "
        "database (BugBuster --sylph_read_seq_id corrects it). Taxonomic abundance is relative to the profiled "
        "community.",
        "Tier rules and the watchlist are editable tables (assets/arg_priority_tiers.tsv, assets/summary_watchlist.tsv).",
    ]
    if ctx.get("watchlist") is not None:
        caveats.append("Watchlist: the taxa in the table below, from " + str(ctx.get("watchlist_path", "")) + ". The "
                       "default list holds the ESKAPE pathogens, the bacteria of the WHO Bacterial Priority Pathogens "
                       "List 2024 and of the CDC Antibiotic Resistance Threats Report 2019, and a few other species of "
                       "clinical interest; each row names its source. A species row matches the species, a genus "
                       "row every species of the genus, GTDB suffixes removed.")
    methods = ("<ul class='caveats'>" + "".join(f"<li>{html.escape(c)}</li>" for c in caveats) + "</ul>"
               + "<h3>Inputs</h3><pre>" + html.escape("\n".join(ctx["inputs"])) + "</pre>")
    if ctx.get("watchlist") is not None:
        methods += "<h3>Watchlist</h3>" + df_to_html(ctx["watchlist"], table_id="t_wldef")
    if ctx["versions"]:
        methods += "<h3>Software versions</h3><pre>" + html.escape(ctx["versions"]) + "</pre>"
    section("methods", "Methods and caveats", methods)

    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title>
<style>{HTML_CSS}</style></head><body>
<header><h1>{html.escape(title)}</h1><div class="muted">BugBuster results summary · generated {date.today().isoformat()}
· {s['n_samples']} samples from {s['n_batches']} results director{'y' if s['n_batches'] == 1 else 'ies'}</div>
<div class="cards">{card_html}</div>{attention}</header>
<nav>{''.join(nav)}</nav><main>{''.join(sections)}</main><script>{HTML_JS}</script></body></html>"""
    path = outdir / "results_summary.html"
    path.write_text(doc)
    return path
