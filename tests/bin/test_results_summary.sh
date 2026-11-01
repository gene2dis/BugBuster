#!/bin/bash
#
# Regression tests for bin/results_summary.py on the fixtures in
# tests/data/summary/ (tests/bin/make_summary_fixtures.sh) and the K-MARVEL
# tables in tests/data/kmarvel/ (sampleA: TEM + CTX-M, sampleB: no ARGs):
#   - read pairs = Reads_report.csv reads / 2; low-depth flag
#   - RPM / RPKM against a hand calculation; tiers from the default rule table
#   - a gene not found in a low-depth sample is NA_low_depth, not 0
#   - sylph ANI interval below 95% -> species_call_uncertain
#   - two results directories merge (taxa outer-joined), a sample in both fails
#   - CARD snps.txt check of variant-model mutations, gene names and tier rules
#     (cfr genes are tier 1)
#   - RGI bwt (tests/data/summary/rgi/): genes at the --rgi-min-coverage cutoff
#     (best reference sequence), reads summed per ARO and gene, k-mer host
#     check against the watchlist, the K-MARVEL vs RGI comparison and overlap
#     counts; an RGI-only directory makes RGI the primary ARG method
#   - depth classes, --taxa-per-sample, heatmap rows by phylum then A-Z
#   - HTML is self-contained (no external URLs)
#   - ARGs-OAP (tests/data/summary/args_oap/): SARG subtypes and types in copies
#     per cell, the cell estimates in samples_qc, low-depth and not-normalised cells
#   - read functional and MicrobeCensus (tests/data/summary/functional/): CPGE
#     tables against the source, COG names, NA_no_ags without MicrobeCensus
#   - README.md / README.docx (well-formed WordprocessingML) / LLM_GUIDE.md, the
#     stages read from an invalid software_versions.yml, "ran, no table found"
#   - contig ARGs on geNomad's plasmids / viruses / proviruses
#     (tests/data/summary/contig/): arg_contig_mobile per assembly and class
#   - every emitted column documented; taxonomy NA_low_depth and prevalence
#     denominator; one tier per gene across methods; CRP-1 keeps its name
#   - stage tables (tests/data/summary/stages/: a per-sample assembly and a
#     co-assembly results tree trimmed from a real run, samples renamed):
#     assembly statistics (BBMap's count/length N/L50), reads on contigs and
#     genes, MAGs with MIMAG tiers (rRNA/tRNA from Bakta, checkm2_only cap),
#     MAG depth (length-weighted / padded bin depth files), contig ARGs with
#     taxon and bin, ARG gene abundance keyed on the assembly (repeated gene
#     ids), bin ARG clusters (header parsing, sharing), contig function tables
#     (CAZy names, cazy_dbcan apart), ARG RPKG, the merge of both trees, and
#     the catalog's processes exist in the pipeline
#
# Python runs inside the image pinned in modules/local/results_summary/main.nf
# (docker; set RS_TEST_HOST=1 to use the host python3 instead).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
SUMMARY="${REPO_DIR}/bin/results_summary.py"
FIX="${REPO_DIR}/tests/data/summary"
KM="${REPO_DIR}/tests/data/kmarvel"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

SUMMARY_IMG=$(grep -o "'community.wave.seqera.io/library/python_pandas[^']*'" \
    "${REPO_DIR}/modules/local/results_summary/main.nf" | tr -d "'" | head -1)
# REPORT_TEST_IMG=<image>: run the scripts in that image instead of the pinned one (to try an image before pinning it)
SUMMARY_IMG=${REPORT_TEST_IMG:-${SUMMARY_IMG}}
py() {
    # python3 in the module's image, with the repo and the scratch dir at their host paths
    if [ "${RS_TEST_HOST:-0}" = "1" ]; then
        PYTHONDONTWRITEBYTECODE=1 python3 "$@"
    else
        docker run --rm -i -u "$(id -u):$(id -g)" -e HOME=/tmp -e MPLCONFIGDIR=/tmp -e PYTHONDONTWRITEBYTECODE=1 \
            -v "${REPO_DIR}:${REPO_DIR}:ro" -v "${TMP_DIR}:${TMP_DIR}" -w "$(pwd)" \
            "${SUMMARY_IMG}" python3 "$@"
    fi
}

PASS=0
FAIL=0
check() {
    local desc="$1"
    shift
    if "$@"; then
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${desc}"
        FAIL=$((FAIL + 1))
    fi
}
fails() {
    ! "$@" > /dev/null 2>&1
}

# make_results <dir> <samples...>: a BugBuster results directory for the samples
make_results() {
    local dir="$1"
    shift
    mkdir -p "${dir}/01_quality_control/summary" "${dir}/02_taxonomy/tables" "${dir}/02_taxonomy/sylph"
    head -1 "${FIX}/Reads_report.csv" > "${dir}/01_quality_control/summary/Reads_report.csv"
    for s in "$@"; do
        grep "^${s}," "${FIX}/Reads_report.csv" >> "${dir}/01_quality_control/summary/Reads_report.csv"
        mkdir -p "${dir}/05_arg_prediction/read_level/kmarvel/${s}"
        cp "${KM}/${s}_kmarvel.tsv" "${dir}/05_arg_prediction/read_level/kmarvel/${s}/"
        [ -e "${FIX}/${s}_sylph_gtdb.tsv" ] && cp "${FIX}/${s}_sylph_gtdb.tsv" "${dir}/02_taxonomy/sylph/"
        if [ -e "${FIX}/rgi/${s}_rgi_bwt.allele_mapping_data.txt" ]; then
            mkdir -p "${dir}/05_arg_prediction/read_level/rgi/${s}"
            cp "${FIX}/rgi/${s}_rgi_bwt.allele_mapping_data.txt" "${dir}/05_arg_prediction/read_level/rgi/${s}/"
        fi
        if [ -e "${FIX}/rgi/${s}_rgi_kmer_61mer_analysis.gene.txt" ]; then
            mkdir -p "${dir}/05_arg_prediction/read_level/rgi_kmer/${s}"
            cp "${FIX}/rgi/${s}_rgi_kmer_61mer_analysis.gene.txt" "${dir}/05_arg_prediction/read_level/rgi_kmer/${s}/"
        fi
    done
    # the tables hold sampleA only (sampleB has no sylph genome, so a run of
    # sampleB alone has no taxonomy tables)
    case " $* " in *" sampleA "*) cp "${FIX}"/sylph_gtdb_*.tsv "${dir}/02_taxonomy/tables/" ;; esac
    # run-level ARGs-OAP and read functional tables, with this directory's samples only
    local oap="${dir}/05_arg_prediction/read_level/args_oap_summary" func="${dir}/07_functional_annotation/summary"
    mkdir -p "${oap}" "${func}"
    for t in args_oap_long args_oap_cells; do
        head -1 "${FIX}/args_oap/${t}.tsv" > "${oap}/${t}.tsv"
        for s in "$@"; do grep "^${s}"$'\t' "${FIX}/args_oap/${t}.tsv" >> "${oap}/${t}.tsv" || true; done
    done
    for t in read_function_abundance read_annotated_fraction read_sample_summary; do
        head -1 "${FIX}/functional/${t}.tsv" > "${func}/${t}.tsv"
        for s in "$@"; do grep "^${s}"$'\t' "${FIX}/functional/${t}.tsv" >> "${func}/${t}.tsv" || true; done
    done
    for s in "$@"; do
        if [ -e "${FIX}/functional/${s}.ags.tsv" ]; then
            mkdir -p "${dir}/07_functional_annotation/microbecensus/${s}"
            cp "${FIX}/functional/${s}.ags.tsv" "${dir}/07_functional_annotation/microbecensus/${s}/"
        fi
    done
    # software_versions.yml as published by older runs: the (removed) KARGVA entry indented with its heredoc
    # terminator left in (tab-indented module) and RGI's multi-line junk, so it is not valid YAML
    mkdir -p "${dir}/pipeline_info"
    printf '       "KARGVA":\n           kargva: 1.0\n       END_VERSIONS\n"KARGA":\n    karga: 1.1\n"KMARVEL":\n    kmarvel: 1.0.0\n"ARGS_OAP":\n    args_oap: 4.0.0\n"RGI_BWT":\n    rgi:\n\n\nunknown\n    kma:  Invalid option:\t--version\n"QC:FASTP":\n    fastp: 1.0.1\n"READ_FUNCTIONAL:WOLTKA_CLASSIFY":\n    woltka: 0.1.7\n"TAXONOMY:SYLPH_PROFILE":\n    sylph: 1.0.0\n' \
        > "${dir}/pipeline_info/software_versions.yml"
}

echo "=== results_summary.py test suite ==="
[ "${RS_TEST_HOST:-0}" = "1" ] && echo "python: host" || echo "image: ${SUMMARY_IMG}"
cd "${TMP_DIR}"
make_results run1/results sampleA sampleB
mkdir -p run1/results/05_arg_prediction/contig_level/summary
cp "${FIX}/contig/Contig_tax_and_arg_prediction.tsv" run1/results/05_arg_prediction/contig_level/summary/

check "single run with metadata, grouping and ordering" \
    py "${SUMMARY}" --results run1/results --outdir out --metadata "${FIX}/metadata.tsv" \
        --group-by subject --order-by visit --min-clean-pairs 10000
for f in samples_qc taxonomy_long taxonomy_species_wide taxonomy_species_prevalence watchlist arg_long \
         arg_gene_rpkm arg_presence arg_drug_class_rpkm arg_priority_panel group_arg_detection; do
    check "table ${f}.tsv written" test -s "out/tables/${f}.tsv"
done
check "Excel workbook written" test -s out/results_summary.xlsx
check "HTML report written" test -s out/results_summary.html
check "HTML has no external URLs" fails grep -Eq '(src|href)="https?://' out/results_summary.html
check "QC, taxonomy and ARG figures written" test -s out/figures/qc_read_fate.png -a -s out/figures/taxonomy_species_heatmap.svg -a -s out/figures/arg_priority_heatmap.pdf

check "table contents" py - "${KM}" <<'PYEOF'
import sys
import pandas as pd
km = sys.argv[1]
qc = pd.read_csv("out/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert qc.at["sampleA", "clean_pairs"] == 60000, qc.at["sampleA", "clean_pairs"]
assert qc.at["sampleB", "raw_pairs"] == 2500
assert not qc.at["sampleA", "low_depth"] and qc.at["sampleB", "low_depth"]
assert qc.at["sampleB", "arg_status"] == "NA_low_depth"
assert qc.at["sampleA", "subject"] == "S1"

arg = pd.read_csv("out/tables/arg_long.tsv", sep="\t")
ref = pd.read_csv(f"{km}/sampleA_kmarvel.tsv", sep="\t")
for name, gene, tier, cat in (("CTX-M-107", "CTX-M", 2, "ESBL"), ("TEM-117", "TEM", 2, "broad_spectrum_bla")):
    row = arg[arg["aro_full_name"] == name].iloc[0]
    assert (row["gene"], row["tier"], row["category"]) == (gene, tier, cat), row
    r = ref[ref["aro_full_name"] == name].iloc[0]
    rpm = r["supporting_reads"] / 60000 * 1e6
    rpkm = rpm / (len(r["reference_sequence"]) * 3 / 1000)
    assert abs(row["rpm"] - rpm) < 1e-6 and abs(row["rpkm"] - rpkm) < 1e-6, (row["rpm"], rpm, row["rpkm"], rpkm)

pres = pd.read_csv("out/tables/arg_presence.tsv", sep="\t", index_col=0)
assert str(pres.at["CTX-M", "sampleA"]) == "1"
assert pres.at["CTX-M", "sampleB"] == "NA_low_depth"

tax = pd.read_csv("out/tables/taxonomy_long.tsv", sep="\t").set_index("taxon_id")
assert bool(tax.at["GCA_047199095.1", "species_call_uncertain"])
assert not bool(tax.at["GCF_000009045.1", "species_call_uncertain"])
assert tax.at["GCF_000009045.1", "est_read_pairs"] == 12507

grp = pd.read_csv("out/tables/group_arg_detection.tsv", sep="\t")
assert set(grp["pattern"]) == {"V1:+ V2:NA"}, grp["pattern"]
PYEOF

check "contig ARGs on plasmids, viruses and proviruses (arg_contig_mobile)" py - <<'PYEOF'
import pandas as pd
mob = pd.read_csv("out/tables/arg_contig_mobile.tsv", sep="\t")
assert list(mob.columns[:3]) == ["assembly", "batch", "arg_class"] and set(mob["batch"]) == {"run1"}, mob
mob = mob.set_index(["assembly", "arg_class"])
assert list(mob.index) == [("sampleA", "tetracycline"), ("sampleA", "bacitracin"), ("sampleA", "multidrug")], mob
tet = mob.loc[("sampleA", "tetracycline")]
assert (tet["arg_genes"], tet["arg_genes_on_plasmid"], tet["arg_genes_on_virus"], tet["arg_genes_in_provirus"]) \
    == (2, 1, 0, 0) and tet["pct_mobile"] == 50.0, tet
assert mob.loc[("sampleA", "bacitracin"), "arg_genes_on_virus"] == 1
assert mob.loc[("sampleA", "multidrug"), "arg_genes_in_provirus"] == 1
md = open("out/README.md").read()
assert "ICE" in md or "integrative conjugative" in md
for stage in ("| Contig taxonomy and ARGs", "| Plasmid / virus classification of contigs (geNomad)"):
    row = next(l for l in md.splitlines() if l.startswith(stage))
    assert "| ran |" in row, row
PYEOF

# two merged co-assembly batches: both contig tables say 'coassembly', so the batch keeps them apart
make_results co1/results sampleA
make_results co2/results sampleB
for b in co1 co2; do
    mkdir -p "${b}/results/05_arg_prediction/contig_level/summary"
    sed 's/^sampleA\t/coassembly\t/' "${FIX}/contig/Contig_tax_and_arg_prediction.tsv" \
        > "${b}/results/05_arg_prediction/contig_level/summary/Contig_tax_and_arg_prediction.tsv"
done
check "merge of two co-assembly batches" \
    py "${SUMMARY}" --results co1/results co2/results --outdir co_out --min-clean-pairs 10000 --no-excel
check "arg_contig_mobile keeps the two co-assemblies apart (batch in the key)" py - <<'PYEOF'
import pandas as pd
mob = pd.read_csv("co_out/tables/arg_contig_mobile.tsv", sep="\t")
tet = mob[mob["arg_class"] == "tetracycline"]
assert sorted(tet["batch"]) == ["co1", "co2"] and set(tet["assembly"]) == {"coassembly"}, tet
assert list(tet["arg_genes"]) == [2, 2], tet
PYEOF

check "ARGs-OAP tables: SARG subtypes and types in copies per cell, cell estimates" py - <<'PYEOF'
import pandas as pd
sl = pd.read_csv("out/tables/arg_sarg_long.tsv", sep="\t")
assert sorted(sl["sarg_subtype"]) == ["multidrug__adeF", "novobiocin__novA", "polymyxin__rosB",
                                      "tetracycline__tetA(48)"] and set(sl["sample"]) == {"sampleA"}, sl
tet = sl[sl["sarg_subtype"] == "tetracycline__tetA(48)"].iloc[0]
# 4 reads (2 from the singletons), 0.548105 copies over 16 cells (pairs + singletons)
assert tet["read_count"] == 4 and abs(tet["copies_per_cell"] - 0.548104956 / 16) < 1e-9
assert {"subject", "visit", "batch"} <= set(sl.columns)
typ = pd.read_csv("out/tables/arg_sarg_type_cpc.tsv", sep="\t", index_col=0)
assert abs(float(typ.at["tetracycline", "sampleA"]) - 0.548104956 / 16) < 1e-9
# sampleB is low depth (--min-clean-pairs 10000) with no SARG hit: unknown, not absent
assert typ.at["novobiocin", "sampleB"] == "NA_low_depth"
sub = pd.read_csv("out/tables/arg_sarg_subtype_cpc.tsv", sep="\t", index_col=[0, 1])
assert ("multidrug", "multidrug__adeF") in sub.index
qc = pd.read_csv("out/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert qc.at["sampleA", "args_oap_ncell"] == 16 and qc.at["sampleA", "args_oap_n16S"] == 32
assert qc.at["sampleA", "sarg_n_subtypes"] == 4
assert abs(qc.at["sampleA", "sarg_copies_per_cell"] - pd.to_numeric(typ["sampleA"]).sum()) < 1e-12
assert abs(qc.at["sampleA", "sarg_copies_per_cell"] - 0.895062 / 16) < 1e-6
PYEOF
check "ARGs-OAP figure written" test -s out/figures/arg_sarg_type_heatmap.svg

check "read functional and MicrobeCensus: CPGE tables, COG names, samples_qc columns" py - "${FIX}" <<'PYEOF'
import sys
import pandas as pd
fix = sys.argv[1]
ko = pd.read_csv("out/tables/function_read_ko.tsv", sep="\t", index_col=0)
src = pd.read_csv(f"{fix}/functional/read_function_abundance.tsv", sep="\t")
row = src[(src["sample_id"] == "sampleA") & (src["ontology"] == "ko") & (src["accession"] == "K00001")].iloc[0]
assert abs(float(ko.at["K00001", "sampleA"]) - row["abundance_cpge"]) < 1e-9
assert ko.at["K00001", "description"].startswith("E1.1.1.1")
# sampleB is low depth and has no MicrobeCensus: no CPGE, absence unknown
assert set(ko["sampleB"]) == {"NA_low_depth"}
cog = pd.read_csv("out/tables/function_read_cog.tsv", sep="\t", index_col=0)
assert cog.at["E", "description"] == "Amino acid transport and metabolism"
long = pd.read_csv("out/tables/function_read_long.tsv", sep="\t")
assert set(long["sample"]) == {"sampleA", "sampleB"} and "abundance_tpm" not in long
assert long.loc[long["ontology"] == "cog", "term_description"].notna().all()
qc = pd.read_csv("out/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert qc.at["sampleA", "mc_genome_equivalents"] == 2.0 and qc.at["sampleA", "mc_avg_genome_size_bp"] == 3e6
assert qc.at["sampleA", "func_reads_assigned"] == 20 and qc.at["sampleA", "func_fraction_annotated"] == 0.9
assert pd.isna(qc.at["sampleB", "mc_genome_equivalents"])
PYEOF
check "COG figure written" test -s out/figures/function_cog_heatmap.png
check "--min-clean-pairs 100: a sample without MicrobeCensus is NA_no_ags in the CPGE tables" \
    py "${SUMMARY}" --results run1/results --outdir lowmin --min-clean-pairs 100 --no-excel
check "NA_no_ags" py - <<'PYEOF'
import pandas as pd
ko = pd.read_csv("lowmin/tables/function_read_ko.tsv", sep="\t", index_col=0)
assert set(ko["sampleB"]) == {"NA_no_ags"} and pd.to_numeric(ko["sampleA"]).gt(0).any()
PYEOF

check "README.md, README.docx and LLM_GUIDE.md" py - <<'PYEOF'
import os, re, zipfile
from xml.dom import minidom
md = open("out/README.md").read()
for text in ("## What ran", "## Tables", "## Figures", "## Units and missing values", "## ARG reporting tiers",
             "`tables/arg_long.tsv`", "`tables/arg_sarg_long.tsv`", "## Software versions", "ARGs-OAP"):
    assert text in md, text
# stage table: the fixture's versions file lists RGI_BWT and the RGI table is there; no assembly
row = next(l for l in md.splitlines() if l.startswith("| Read-level ARGs: RGI"))
assert "| ran |" in row, row
assert re.search(r"\| Assembly \(MEGAHIT\) \| not run \|", md)
# versions read line by line from the invalid YAML
assert "| KARGVA | kargva | 1.0 |" in md and "| RGI_BWT | kma | unknown |" in md
z = zipfile.ZipFile("out/README.docx")
assert {"[Content_Types].xml", "word/document.xml", "word/styles.xml", "word/numbering.xml"} <= set(z.namelist())
for name in z.namelist():
    minidom.parseString(z.read(name))
doc = z.read("word/document.xml").decode()
assert 'w:val="Heading2"' in doc and "What ran" in doc and "<w:tbl>" in doc
guide = open("out/LLM_GUIDE.md").read()
for text in ("## Run facts", "n_samples: 2", "## Pitfalls", "### arg_long.tsv", "### arg_sarg_long.tsv", "### samples_qc.tsv",
             "`tables/samples_qc.tsv`", "## Figures (figures/)", "### Binning and MAGs", "## Recipes"):
    assert text in guide, text
# the summary folder is what gets shared: paths relative to it, no local path, inputs by file name
html = open("out/results_summary.html").read()
for doc in (md, guide, html):
    assert os.getcwd() not in doc and "00_summary/" not in doc and "results_dir" not in doc
    assert not re.search(r"--[a-z-]+ /\w", doc)  # a flag followed by an absolute path
assert "--arg-priority arg_priority_tiers.tsv" in md and "--results: " in md
assert "05_arg_prediction" not in guide
# every column row of the guide has a meaning
cols = re.findall(r"^\| `([^`]+)` \| (.*) \|$", guide, flags=re.M)
assert cols and all(meaning.strip() for _, meaning in cols), [c for c, m in cols if not m.strip()][:5]
PYEOF
check "Excel README sheet lists every table written" py - <<'PYEOF'
import os
import pandas as pd
readme = pd.read_excel("out/results_summary.xlsx", sheet_name="README")
tables = {f[:-4][:31] for f in os.listdir("out/tables")}
assert set(readme["table"]) == tables, set(readme["table"]) ^ tables
PYEOF

check "RGI genes at the coverage cutoff, k-mer hosts and the method comparison" py - <<'PYEOF'
import pandas as pd
arg = pd.read_csv("out/tables/arg_long.tsv", sep="\t")
rgi = arg[arg["method"] == "rgi"].set_index("aro_full_name")
# KPC-3 (30%) and CTX-M-15 (50%) are below the 80% cutoff
assert set(rgi.index) == {"KPC-2", "TEM-1", "poxtA", "cfrE"}, set(rgi.index)
kpc = rgi.loc["KPC-2"]
# two reference sequences: best coverage 100%, reads summed (30 + 10); read pairs = reads / 2
assert kpc["rgi_mapped_reads"] == 40 and kpc["assembly_coverage"] == 100 and kpc["supporting_reads"] == 20
rpkm = 20 / 60000 * 1e6 / 0.882
assert abs(kpc["rpkm"] - rpkm) < 1e-6
# the gene sums all its AROs, KPC-3's 4 reads (2 pairs) included
assert abs(kpc["gene_rpkm"] - (rpkm + 2 / 60000 * 1e6 / 0.882)) < 1e-6
assert kpc["rgi_kmer_hosts"].startswith("Klebsiella pneumoniae (25)") and kpc["rgi_host_check"] == "watchlist_host"
assert rgi.at["poxtA", "rgi_host_check"] == "non_watchlist_host"
assert rgi.at["poxtA", "rgi_best_ref"] == "WildCARD: Blautia producta"
assert rgi.at["cfrE", "rgi_host_check"] == "no_kmer_host" and rgi.at["cfrE", "tier"] == 1
assert rgi.at["TEM-1", "rgi_host_check"] == "watchlist_host"
comp = pd.read_csv("out/tables/arg_method_comparison.tsv", sep="\t").set_index("gene")
assert comp.at["TEM", "detected_by"] == "K-MARVEL + RGI"
assert comp.at["CTX-M", "detected_by"] == "K-MARVEL only" and bool(comp.at["CTX-M", "rgi_below_cutoff"])
assert comp.at["KPC", "detected_by"] == "RGI only" and comp.at["KPC", "rgi_read_pairs"] == 22
ov = pd.read_csv("out/tables/arg_method_overlap.tsv", sep="\t")
row = ov[(ov["unit"] == "genes") & (ov["level"] == "all")].iloc[0]
got = (row["K-MARVEL + RGI"], row["K-MARVEL only"], row["RGI only"], row["K-MARVEL only, RGI below cutoff"], row["total"])
assert got == (1, 1, 3, 1, 5), got
qc = pd.read_csv("out/tables/samples_qc.tsv", sep="\t").set_index("sample")
got = (qc.at["sampleA", "n_genes_kmarvel"], qc.at["sampleA", "n_genes_rgi"], qc.at["sampleA", "n_genes_all_methods"])
assert got == (2, 4, 1), got
assert qc.at["sampleB", "n_genes_rgi"] == 0
# depth classes: --min-clean-pairs 10000 and the default marks (10M, 20M)
assert qc.at["sampleA", "depth_class"] == "10k-10M" and qc.at["sampleB", "depth_class"] == "<10k"
html = open("out/results_summary.html").read()
assert "Comparison of ARG methods" in html and "possibly a related gene of a commensal" in html
assert "samples with &gt;=10M clean pairs" in html
# the reading guide comes before the first figure, and explains tiers and both methods
assert html.index("Reading this report</h2>") < html.index("<figure>")
assert "ARG reporting tiers" in html and "RGI k-mer (pathogen of origin)" in html
import os
for gone in ("taxonomy_diversity", "taxonomy_pcoa", "group_timeline", "arg_tier1_timeline"):
    assert not os.path.exists(f"out/figures/{gone}.png"), gone
PYEOF
check "method comparison figures written" test -s out/figures/arg_method_upset.svg -a -s out/figures/arg_method_by_category.png
documented() {
    # documented <outdir> <args...>: a run whose every emitted column is in summary_docs.TABLE_DOCS
    local out="$1"
    shift
    py "${SUMMARY}" --outdir "${out}" --no-excel "$@" > "${out}.log" 2>&1 || return 1
    ! grep -q "columns without documentation" "${out}.log"
}
check "every column documented (metadata, grouping, every method)" \
    documented docs1 --results run1/results --metadata "${FIX}/metadata.tsv" --group-by subject --order-by visit \
        --min-clean-pairs 10000
check "no stale column documentation (only columns this fixture cannot produce are unused)" py - "${REPO_DIR}/bin" <<'PYEOF'
import os, sys
import pandas as pd
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
from summary_docs import unused_doc_columns
tables = {}
for f in sorted(os.listdir("docs1/tables")):
    name = f[:-4]
    df = pd.read_csv(f"docs1/tables/{f}", sep="\t")
    tables[name] = (df, False)
unused = set(unused_doc_columns(tables, "subject"))
# --card-prevalence columns and mean_read_length (fastp JSON) have no fixture
allowed = {("arg_long", "card_plasmid_associated"), ("arg_long", "candidate_hosts"), ("samples_qc", "mean_read_length"),
           # the assembly and contig functional columns (tested on tests/data/summary/stages below)
           ("samples_qc", "pct_reads_on_genes"), ("samples_qc", "func_contig_genes_annotated_pct"),
           ("samples_qc", "func_contig_annotated_by_abundance_pct")}
assert unused <= allowed, sorted(unused - allowed)
PYEOF

# metadata columns named like common study variables are kept as they are
printf 'sample\tgroup\ttype\tclass\tdescription\nsampleA\tg1\tt1\tc1\td1\nsampleB\tg2\tt2\tc2\td2\n' > common_meta.tsv
check "metadata with group / type / class / description columns" \
    documented commonmeta --results run1/results --metadata common_meta.tsv --group-by group --min-clean-pairs 10000
check "those columns are not renamed" py - <<'PYEOF'
import pandas as pd
qc = pd.read_csv("commonmeta/tables/samples_qc.tsv", sep="\t")
assert {"group", "type", "class", "description"} <= set(qc.columns), qc.columns
al = pd.read_csv("commonmeta/tables/arg_long.tsv", sep="\t")
assert {"group", "gene"} <= set(al.columns), al.columns
PYEOF

# low depth in the taxonomy tables: sampleA below --min-clean-pairs
check "--min-clean-pairs 100000 makes sampleA low depth" \
    py "${SUMMARY}" --results run1/results --outdir lowA --min-clean-pairs 100000 --no-excel
check "taxonomy: a taxon not found in a low-depth sample is NA_low_depth; prevalence leaves it out" py - <<'PYEOF'
import pandas as pd
wl = pd.read_csv("lowA/tables/watchlist.tsv", sep="\t").set_index("taxon")
assert wl.at["Klebsiella pneumoniae", "sampleA"] == "NA_low_depth", wl.loc["Klebsiella pneumoniae"]
prev = pd.read_csv("lowA/tables/taxonomy_species_prevalence.tsv", sep="\t")
assert (prev["n_samples_evaluable"] == prev["n_samples_detected"]).all()
assert (prev["pct_samples_detected"] == 100).all()
PYEOF

check "--rgi-min-coverage 20 detects the low-coverage genes" \
    py "${SUMMARY}" --results run1/results --outdir cut20 --min-clean-pairs 10000 --rgi-min-coverage 20 --no-excel
check "cutoff 20: KPC-3 and CTX-M-15 are RGI calls" py - <<'PYEOF'
import pandas as pd
arg = pd.read_csv("cut20/tables/arg_long.tsv", sep="\t")
assert {"KPC-3", "CTX-M-15"} <= set(arg.loc[arg["method"] == "rgi", "aro_full_name"])
comp = pd.read_csv("cut20/tables/arg_method_comparison.tsv", sep="\t").set_index("gene")
assert comp.at["CTX-M", "detected_by"] == "K-MARVEL + RGI"
PYEOF
check "--taxa-per-sample 1 draws one species row" \
    py "${SUMMARY}" --results run1/results --outdir top1 --min-clean-pairs 10000 --taxa-per-sample 1 --no-excel
check "species heatmap caption gives one row" grep -q "most abundant species taxa (1 row)" top1/results_summary.html

# Edge cases that used to crash the report: no tier 1-2 gene in any group, an empty watchlist,
# metadata columns named like summary columns
printf 'match\tpattern\ttier\tcategory\tnote\nname\t.\t3\tall_tier3\ttest: every gene tier 3\n' > tier3.tsv
printf 'taxon\trank\tnote\n' > empty_watchlist.tsv
printf 'sample\tsubject\tvisit\tclean_pairs\ttier\nsampleA\tS1\tV1\t5\tx\nsampleB\tS1\tV2\t6\ty\n' > clash_meta.tsv
check "no tier 1-2 gene with --group-by, empty watchlist, clashing metadata columns" \
    py "${SUMMARY}" --results run1/results --outdir edge --min-clean-pairs 10000 --no-excel \
        --metadata clash_meta.tsv --group-by subject --order-by visit --arg-priority tier3.tsv \
        --watchlist empty_watchlist.tsv
check "edge cases: metadata renamed, computed columns kept, empty tables typed" py - <<'PYEOF'
import pandas as pd
qc = pd.read_csv("edge/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert qc.at["sampleA", "clean_pairs"] == 60000 and str(qc.at["sampleA", "meta_clean_pairs"]) == "5"
assert qc.at["sampleA", "meta_tier"] == "x"
grp = pd.read_csv("edge/tables/group_arg_detection.tsv", sep="\t")
assert grp.empty and "tier" in grp.columns
wl = pd.read_csv("edge/tables/watchlist.tsv", sep="\t")
assert wl.empty and "n_samples_detected" in wl.columns
PYEOF

# RGI without K-MARVEL: RGI is the primary ARG method, nothing to compare
make_results rgionly/results sampleA sampleB
rm -r rgionly/results/05_arg_prediction/read_level/kmarvel
check "RGI-only results directory" py "${SUMMARY}" --results rgionly/results --outdir rgi1 --min-clean-pairs 10000 --no-excel
check "RGI-only: RGI drives the ARG matrices, no comparison" py - <<'PYEOF'
import os
import pandas as pd
assert set(pd.read_csv("rgi1/tables/arg_long.tsv", sep="\t")["method"]) == {"rgi"}
panel = pd.read_csv("rgi1/tables/arg_priority_panel.tsv", sep="\t", index_col=0)
assert {"KPC", "cfrE", "poxtA"} <= set(panel.index), panel.index
assert not os.path.exists("rgi1/tables/arg_method_comparison.tsv")
assert "Antimicrobial resistance genes (RGI)" in open("rgi1/results_summary.html").read()
PYEOF
check "--arg-primary kmarvel without K-MARVEL output fails" \
    fails py "${SUMMARY}" --results rgionly/results --outdir bad --arg-primary kmarvel --no-excel

# Two results directories; sampleB's has no taxonomy tables
make_results batchA/results sampleA
make_results batchB/results sampleB
check "two results directories merge" \
    py "${SUMMARY}" --results batchA/results batchB/results --outdir merged --min-clean-pairs 10000 --no-excel
check "merged tables match the single run" py - <<'PYEOF'
import pandas as pd
for name in ("arg_long", "arg_presence", "taxonomy_species_wide"):
    a = pd.read_csv(f"out/tables/{name}.tsv", sep="\t")
    b = pd.read_csv(f"merged/tables/{name}.tsv", sep="\t")
    cols = [c for c in a.columns if c not in ("batch", "subject", "visit")]
    pd.testing.assert_frame_equal(a[cols].reset_index(drop=True), b[cols].reset_index(drop=True))
qc = pd.read_csv("merged/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert qc.at["sampleA", "batch"] == "batchA" and qc.at["sampleB", "batch"] == "batchB"
PYEOF
check "--batch-labels names the batches" \
    py "${SUMMARY}" --results batchA/results batchB/results --batch-labels b1 b2 --outdir labelled \
        --min-clean-pairs 10000 --no-excel
check "batch column holds the given labels" py - <<'PYEOF'
import pandas as pd
qc = pd.read_csv("labelled/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert qc.at["sampleA", "batch"] == "b1" and qc.at["sampleB", "batch"] == "b2"
PYEOF
check "--batch-labels with fewer labels than directories fails" \
    fails py "${SUMMARY}" --results batchA/results batchB/results --batch-labels b1 --outdir bad --no-excel
check "duplicate --batch-labels fail" \
    fails py "${SUMMARY}" --results batchA/results batchB/results --batch-labels b1 b1 --outdir bad --no-excel
check "K-MARVEL ran but found no ARG in any sample" \
    py "${SUMMARY}" --results batchB/results --outdir noarg --min-clean-pairs 100 --no-excel
check "no-ARG run: zero calls reported, no ARG matrices" py - <<'PYEOF'
import os
import pandas as pd
qc = pd.read_csv("noarg/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert qc.at["sampleB", "arg_status"] == "ok" and qc.at["sampleB", "n_arg_calls_tier1"] == 0
assert not os.path.exists("noarg/tables/arg_long.tsv")
PYEOF
check "a sample in two results directories fails" \
    fails py "${SUMMARY}" --results run1/results batchA/results --outdir dup --no-excel
check "an unknown --group-by column fails" \
    fails py "${SUMMARY}" --results run1/results --outdir bad --metadata "${FIX}/metadata.tsv" --group-by patient --no-excel

# QC only: taxonomy and ARG sections are skipped
mkdir -p qconly/results/01_quality_control/summary
cp "${FIX}/Reads_report.csv" qconly/results/01_quality_control/summary/
check "QC-only results directory" py "${SUMMARY}" --results qconly/results --outdir qc --no-excel
check "QC-only: only samples_qc.tsv" test "$(ls qc/tables)" = "samples_qc.tsv"

# kraken2's Reads_report.csv holds Classified as a percentage
mkdir -p kraken/results/01_quality_control/summary
printf 'Id,Raw reads,Fastp,Clean reads (contaminants removed),Kraken DB,Unclassified,Classified\nk1,2000,1900,1800,standard-8,14.5,85.5\nk2,2000,1900,1800,standard-8,99.5,0.5\n' \
    > kraken/results/01_quality_control/summary/Reads_report.csv
check "kraken2 Reads_report.csv" py "${SUMMARY}" --results kraken/results --outdir kr --no-excel
check "kraken2 classified percentage becomes a fraction" py - <<'PYEOF'
import pandas as pd
qc = pd.read_csv("kr/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert abs(qc.at["k1", "classified_fraction"] - 0.855) < 1e-9 and abs(qc.at["k2", "classified_fraction"] - 0.005) < 1e-9
PYEOF

check "CARD snps.txt check, gene names, tier rules, heatmap row order and k-mer parsing" py - "${REPO_DIR}/bin" <<'PYEOF'
import sys, tempfile, os
import pandas as pd
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
import summary_args as rs
from summary_common import strip_gtdb
from summary_figures import lineage_order
d = tempfile.mkdtemp()
with open(os.path.join(d, "snps.txt"), "w") as fh:
    fh.write("Accession\tName\tModel Type\tParameter Type\tMutations\tCARD Short Name\tsource\tcitation\n")
    fh.write("3003308\tparC\tprotein variant model\tsingle resistance variant\tS80I\tx\tx\tx\n")
    fh.write("3003381\tsoxR\tprotein variant model\tmultiple resistance variants\tT38S,G74R\tx\tx\tx\n")
    fh.write("3003378\tmarR\tprotein overexpression model\tnonsense mutation\tE31Ter\tx\tx\tx\n")
snps = rs.load_snps(d)
assert rs.snp_status("ARO:3003308", "A10T;S80I", snps) == ("card_snp_present", "S80I")
assert rs.snp_status("ARO:3003308", "A10T;S80R", snps) == ("no_card_snp", "")
assert rs.snp_status("ARO:3003381", "T38S", snps) == ("no_card_snp", "")
assert rs.snp_status("ARO:3003381", "G74R;T38S", snps)[0] == "card_snp_present"
assert rs.snp_status("ARO:3003378", "E31*", snps)[0] == "card_snp_present"
assert rs.snp_status("ARO:3000001", "A1B", snps) == ("no_snp_definition", "")
assert rs.gene_name("KPC-2", "KPC beta-lactamase") == "KPC"
assert rs.gene_name("OXA-181", "OXA-48-like beta-lactamase") == "OXA-48-like"
assert rs.gene_name("MCR-9.1", "MCR phosphoethanolamine transferase") == "MCR-9"
assert rs.gene_name("QnrB19", "quinolone resistance protein (qnr)") == "QnrB"
assert rs.gene_name("tet(Q)", "tetracycline-resistant ribosomal protection protein") == "tet(Q)"
assert strip_gtdb("Enterococcus_B faecium") == "Enterococcus faecium"
rules = rs.load_rules(os.path.join(sys.argv[1], "..", "assets", "arg_priority_tiers.tsv"))
def tier(name, family, mech="antibiotic target alteration"):
    row = {"aro_full_name": name, "gene_family": family, "aro_id": "ARO:0", "model": "protein homolog model",
           "resistance_mechanism": mech}
    return rs.assign_tier(row, rules, None)[0]
# a name rule wins over every later family rule: cfr genes must reach the tier 1 Cfr family rule
assert tier("cfrE", "Cfr 23S ribosomal RNA methyltransferase") == 1
assert tier("cfr(B)", "Cfr 23S ribosomal RNA methyltransferase") == 1
assert tier("qacEdelta1", "major facilitator superfamily (MFS) antibiotic efflux pump", "antibiotic efflux") == 2
assert tier("poxtA", "Miscellaneous ABC-F subfamily ATP-binding cassette ribosomal protection proteins") == 1
# vanG (common on the chromosome of gut Clostridia) is tier 3; the other Van ligases stay tier 1
assert tier("vanG", "glycopeptide resistance gene cluster;Van ligase") == 3
for ligase in ("vanA", "vanB", "vanD"):
    assert tier(ligase, "glycopeptide resistance gene cluster;Van ligase") == 1, ligase
# chromosomal genes of Actinomycetota are tier 3: RbpA, and the vanO cluster (Rhodococcus); vanR of other
# clusters stays an accessory gene (tier 2)
assert tier("RbpA", "RbpA bacterial RNA polymerase-binding protein", "antibiotic target protection") == 3
assert tier("vanO", "Van ligase;glycopeptide resistance gene cluster") == 3
assert tier("vanR gene in vanO cluster", "glycopeptide resistance gene cluster;vanR") == 3
assert tier("vanR gene in vanA cluster", "glycopeptide resistance gene cluster;vanR") == 2
assert rs.gene_name("vanR gene in vanO cluster", "glycopeptide resistance gene cluster;vanR") \
    == "vanR gene in vanO cluster"
# heatmap rows: grouped by phylum (largest first), alphabetical within each phylum, GTDB suffixes ignored
tax = pd.DataFrame({"Phylum": ["Bacillota", "Bacillota", "Bacillota", "Bacteroidota"],
                    "Species": ["Zeta b", "alpha_A c", "Beta d", "Aaa e"]}, index=["g1", "g2", "g3", "g4"])
order, phyla, label = lineage_order(["Zeta b", "Aaa e", "Beta d", "alpha_A c"], {"tax": tax, "ranks": ["Phylum", "Species"]},
                                    "Species")
assert order == ["alpha_A c", "Beta d", "Zeta b", "Aaa e"], order
assert phyla == ["Bacillota"] * 3 + ["Bacteroidota"] and label == "Phylum"
# RGI host check: a tie goes to the host rgi_kmer_hosts lists first (alphabetical)
is_watch = rs.watchlist_matcher(pd.DataFrame({"taxon": ["Klebsiella pneumoniae"], "rank": ["species"]}))
assert rs.host_check({"Klebsiella pneumoniae": 5, "Blautia producta": 5}, is_watch, True) == "non_watchlist_host"
assert rs.host_check({"Klebsiella pneumoniae": 6, "Blautia producta": 5}, is_watch, True) == "watchlist_host"
# RGI k-mer prediction strings
assert rs.parse_kmer_prediction("Klebsiella pneumoniae (chromosome or plasmid): 83; Escherichia coli (plasmid): 2; ") \
    == {"Klebsiella pneumoniae": 83, "Escherichia coli": 2}
# a collapsed allele name that is another CARD gene keeps its own name: CRP-1 (CRP beta-lactamase) is not CRP
fams = {"CRP": "resistance-nodulation-cell division (RND) antibiotic efflux pump", "CfxA": "CfxA beta-lactamase"}
assert rs.gene_name("CRP-1", "CRP beta-lactamase", fams) == "CRP-1"
assert rs.gene_name("CfxA2", "CfxA beta-lactamase", fams) == "CfxA"
assert rs.gene_name("CRP-1", "CRP beta-lactamase") == "CRP"
assert rs.norm_list("a; b ;c") == "a;b;c"
# software_versions.yml read line by line: a tab-indented heredoc's terminator and RGI's junk do not break it
from summary_common import parse_versions
v = parse_versions('       "KARGVA":\n           kargva: 1.0\n       END_VERSIONS\n"QC:FASTP":\n    fastp: 1.0.1\n'
                   '"RGI_BWT":\n    rgi:\n\n\nunknown\n    kma:  Invalid option:\t--version\n')
assert v == {"KARGVA": {"kargva": "1.0"}, "FASTP": {"fastp": "1.0.1"},
             "RGI_BWT": {"rgi": "unknown", "kma": "unknown"}}, v
# the catalog names processes as parse_versions keys them: the nf-core aliases (TAXONOMY:KRAKEN2 -> KRAKEN2)
from summary_catalog import stage_status
kv = parse_versions('"TAXONOMY:KRAKEN2":\n    kraken2: 2.1.6\n"TAXONOMY:BRACKEN":\n    bracken: 3.0.1\n')
assert dict((st["key"], s) for st, s in stage_status(kv, set()))["taxonomy"] == "ran, no table found"
# one tier per gene across methods: the first method's tier wins
k = pd.DataFrame({"gene": ["g1"], "tier": [2], "category": ["a"], "tier_note": [""],
                  "resistance_snp_status": [""], "resistance_snps_found": [""]})
r = pd.DataFrame({"gene": ["g1"], "tier": [3], "category": ["b"], "tier_note": [""], "resistance_snp_status": [""]})
assert rs.unify_gene_tiers({"kmarvel": k, "rgi": r}, ["kmarvel", "rgi"]) == ["g1"]
assert r.at[0, "tier"] == 2 and r.at[0, "category"] == "a"
# a K-MARVEL row raised to tier 4 by the reference method gets its CARD SNP check
k4 = pd.DataFrame({"gene": ["g2"], "tier": [2], "category": ["a"], "tier_note": [""], "aro_id": ["ARO:3000001"],
                   "mutation_data": ["S83L"], "resistance_snp_status": [""], "resistance_snps_found": [""]})
r4 = pd.DataFrame({"gene": ["g2"], "tier": [4], "category": ["mutation_dependent"], "tier_note": [""],
                   "resistance_snp_status": ["not_checked_rgi"]})
rs.unify_gene_tiers({"rgi": r4, "kmarvel": k4}, ["rgi", "kmarvel"], {"3000001": [{"S83L"}]})
assert k4.at[0, "tier"] == 4 and k4.at[0, "resistance_snp_status"] == "card_snp_present", k4.iloc[0].to_dict()
PYEOF

# a stage that ran (software_versions.yml) without its table in the results directory is reported
make_results notable/results sampleA
rm -r notable/results/05_arg_prediction/read_level/rgi notable/results/05_arg_prediction/read_level/rgi_kmer
check "RGI_BWT in software_versions.yml without RGI tables" \
    py "${SUMMARY}" --results notable/results --outdir notab --min-clean-pairs 10000 --no-excel
check "stage reported as 'ran, no table found'" \
    grep -q "| Read-level ARGs: RGI bwt + k-mer pathogen of origin (CARD + WildCARD) | ran, no table found |" notab/README.md

# with a profiler, TAXONOMY_REPORT publishes its copy (with the classified fraction) under its own name; the
# summary reads it first, so the QC copy without those columns does not lose the classified fraction
make_results taxrep/results sampleA sampleB
cut -d, -f1-4 "${FIX}/Reads_report.csv" > taxrep/results/01_quality_control/summary/Reads_report.csv
cp "${FIX}/Reads_report.csv" taxrep/results/02_taxonomy/tables/Reads_report_taxonomy.csv
check "02_taxonomy/tables/Reads_report_taxonomy.csv next to a QC copy without Classified" \
    py "${SUMMARY}" --results taxrep/results --outdir taxrep_out --min-clean-pairs 10000 --no-excel
check "classified fraction read from the taxonomy copy" py - <<'PYEOF'
import pandas as pd
qc = pd.read_csv("taxrep_out/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert qc.at["sampleA", "classified_fraction"] == 1.0 and qc.at["sampleB", "classified_fraction"] == 0.0
PYEOF

# no Reads_report.csv: no RPKM, so a gene that was called is blank in the matrices, never 0 (= not detected)
make_results noreads/results sampleA sampleB
rm noreads/results/01_quality_control/summary/Reads_report.csv
check "results directory without Reads_report.csv" \
    py "${SUMMARY}" --results noreads/results --outdir noreads_out --no-excel
check "a call without RPKM is blank in arg_gene_rpkm and 1 in arg_presence" py - <<'PYEOF'
import pandas as pd
arg = pd.read_csv("noreads_out/tables/arg_long.tsv", sep="\t")
called = arg[(arg["method"] == "kmarvel") & arg["gene_representative"]]
rpkm = pd.read_csv("noreads_out/tables/arg_gene_rpkm.tsv", sep="\t", index_col=0, dtype=str, keep_default_na=False)
pres = pd.read_csv("noreads_out/tables/arg_presence.tsv", sep="\t", index_col=0, dtype=str, keep_default_na=False)
assert not called.empty
for g, s in zip(called["gene"], called["sample"]):
    assert rpkm.at[g, s] == "", (g, s, rpkm.at[g, s])
    assert pres.at[g, s] == "1", (g, s, pres.at[g, s])
PYEOF

# tier rules: a tier outside 1-4 fails clearly; the note column is optional
py - "${REPO_DIR}/assets/arg_priority_tiers.tsv" <<'PYEOF'
import sys
import pandas as pd
rules = pd.read_csv(sys.argv[1], sep="\t", comment="#", dtype=str, keep_default_na=False)
bad = pd.concat([rules, pd.DataFrame([{"match": "name", "pattern": "^TEM", "tier": "5", "category": "x", "note": ""}])])
bad.to_csv("tiers_bad.tsv", sep="\t", index=False)
rules.drop(columns="note").to_csv("tiers_nonote.tsv", sep="\t", index=False)
PYEOF
check "a tier 5 rule fails" \
    fails py "${SUMMARY}" --results run1/results --outdir tier5 --arg-priority tiers_bad.tsv --no-excel
py "${SUMMARY}" --results run1/results --outdir tier5 --arg-priority tiers_bad.tsv --no-excel > tier5.log 2>&1 || true
check "the tier 5 error names the tier" grep -q "tier must be 1, 2, 3 or 4" tier5.log
check "a rule table without a note column" \
    py "${SUMMARY}" --results run1/results --outdir nonote --arg-priority tiers_nonote.tsv --no-excel
check "its tier table still lists the rule categories" grep -q "intrinsic_glycopeptide" nonote/README.md

# taxa without a species name are pooled as 'unassigned species': not counted as a species
make_results unas/results sampleA
awk -F'\t' 'BEGIN{OFS="\t"} $1=="GCF_000009045.1"{$NF=""} {print}' "${FIX}/sylph_gtdb_tax_table.tsv" \
    > unas/results/02_taxonomy/tables/sylph_gtdb_tax_table.tsv
check "a taxon without a species name" \
    py "${SUMMARY}" --results unas/results --outdir unas_out --min-clean-pairs 10000 --no-excel
check "n_species and Shannon leave 'unassigned species' out" py - <<'PYEOF'
import numpy as np
import pandas as pd
full = pd.read_csv("out/tables/samples_qc.tsv", sep="\t").set_index("sample")
qc = pd.read_csv("unas_out/tables/samples_qc.tsv", sep="\t").set_index("sample")
wide = pd.read_csv("unas_out/tables/taxonomy_species_wide.tsv", sep="\t", index_col=0)
assert "unassigned species" in wide.index
assert qc.at["sampleA", "n_species"] == full.at["sampleA", "n_species"] - 1
named = wide.loc[[i for i in wide.index if i != "unassigned species"], "sampleA"].astype(float)
p = named[named > 0] / named.sum()
assert abs(qc.at["sampleA", "shannon"] - float(-(p * np.log(p)).sum())) < 1e-9
PYEOF

# ARGs-OAP had no 16S or cell estimate for a sample (stage two stopped): its SARG cells are blank, not 0
make_results nonorm/results sampleA sampleB
sed -i 's/^sampleB\t4000\t1.0\t0.5\t0.0\tno_sarg_hit$/sampleB\t4000\t0.0\t0.0\t0.0\tnot_normalised/' \
    nonorm/results/05_arg_prediction/read_level/args_oap_summary/args_oap_cells.tsv
check "ARGs-OAP sample without a cell estimate" \
    py "${SUMMARY}" --results nonorm/results --outdir nonorm_out --min-clean-pairs 100 --no-excel
check "not normalised: blank SARG cells and counts, cell estimates kept" py - <<'PYEOF'
import pandas as pd
typ = pd.read_csv("nonorm_out/tables/arg_sarg_type_cpc.tsv", sep="\t", index_col=0)
assert typ["sampleB"].isna().all() and (pd.to_numeric(typ["sampleA"]) > 0).all(), typ
qc = pd.read_csv("nonorm_out/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert pd.isna(qc.at["sampleB", "sarg_n_subtypes"]) and qc.at["sampleB", "args_oap_ncell"] == 0
assert "No 16S or cell estimate, so nothing normalised: sampleB" in open("nonorm_out/results_summary.html").read()
PYEOF
# the same sample at low depth (--min-clean-pairs 10000): still blank, not 0 or NA_low_depth
check "low-depth ARGs-OAP sample without a cell estimate" \
    py "${SUMMARY}" --results nonorm/results --outdir nonorm_low --min-clean-pairs 10000 --no-excel
check "not normalised at low depth: blank SARG cells and counts" py - <<'PYEOF'
import pandas as pd
typ = pd.read_csv("nonorm_low/tables/arg_sarg_type_cpc.tsv", sep="\t", index_col=0)
assert typ["sampleB"].isna().all(), typ
qc = pd.read_csv("nonorm_low/tables/samples_qc.tsv", sep="\t").set_index("sample")
assert qc.at["sampleB", "low_depth"] and pd.isna(qc.at["sampleB", "sarg_copies_per_cell"])
assert pd.isna(qc.at["sampleB", "sarg_n_subtypes"])
assert "No 16S or cell estimate, so nothing normalised: sampleB" in open("nonorm_low/results_summary.html").read()
PYEOF

# a 0 cutoff is recorded in the README's inputs, not dropped as if it were unset
check "--rgi-min-coverage 0 --min-clean-pairs 0" \
    py "${SUMMARY}" --results run1/results --outdir zero_out --rgi-min-coverage 0 --min-clean-pairs 0 --no-excel
check "zero settings listed in the README inputs" \
    sh -c 'grep -q "^--rgi-min-coverage 0" zero_out/README.md && grep -q "^--min-clean-pairs 0" zero_out/README.md'

# merged batches where RGI ran on one only: a tier 1 call in the K-MARVEL-only batch is in the header list
make_results kmonly/results sampleB
rm -r kmonly/results/05_arg_prediction/read_level/rgi kmonly/results/05_arg_prediction/read_level/rgi_kmer 2>/dev/null || true
# sampleB's K-MARVEL calls: sampleA's with TEM-117 replaced by KPC-2 (tier 1)
sed 's/^TEM-117\tARO:3000980\tTEM beta-lactamase/KPC-2\tARO:3002312\tKPC beta-lactamase/' "${KM}/sampleA_kmarvel.tsv" \
    > kmonly/results/05_arg_prediction/read_level/kmarvel/sampleB/sampleB_kmarvel.tsv
make_results withrgi/results sampleA
check "merge of a K-MARVEL + RGI batch and a K-MARVEL-only batch" \
    py "${SUMMARY}" --results withrgi/results kmonly/results --outdir mixed_out --min-clean-pairs 10000 --no-excel
check "header lists tier 1 calls of samples outside the method comparison" py - <<'PYEOF'
import pandas as pd
arg = pd.read_csv("mixed_out/tables/arg_long.tsv", sep="\t")
t1 = arg[(arg["tier"] == 1) & (arg["sample"] == "sampleB") & arg["gene_representative"]]
html = open("mixed_out/results_summary.html").read()
assert set(t1["gene"]) == {"KPC"}, set(t1["gene"])
# KPC: RGI only in sampleA (both methods ran), K-MARVEL in sampleB (RGI did not run there)
assert "KPC: 2 samples (sampleA, sampleB)" in html and "K-MARVEL (not every method ran) / RGI only" in html
PYEOF

# ---- stage tables: assembly, MAGs, contig ARGs, bin ARG clusters, contig functions (tests/data/summary/stages)
STG="${FIX}/stages"
make_results st_asm/results sampleA sampleB
cp -r "${STG}/assembly/." st_asm/results/
cp -r "${STG}/coassembly" st_co
check "per-sample assembly tree: every column documented" \
    documented st_asm_out --results st_asm/results --min-clean-pairs 10000
check "co-assembly tree: every column documented" documented st_co_out --results st_co --min-clean-pairs 10000
check "both trees merged (disjoint samples)" \
    documented st_merge --results st_asm/results st_co --batch-labels asm co --min-clean-pairs 10000
check "stage tables against their stage files" py - "${STG}" <<'PYEOF'
import sys
import pandas as pd
stg = sys.argv[1]
t = lambda out, name: pd.read_csv(f"{out}/tables/{name}.tsv", sep="\t", keep_default_na=False)  # noqa: E731
# assembly statistics: BBMap's 'contig N/L50: 1683/8.588 Kbp' is count / length
a = t("st_asm_out", "assembly_stats").set_index("assembly")
assert a.at["sampleA", "l50"] == 1683 and a.at["sampleA", "n50_kb"] == 8.588, a.loc["sampleA"]
assert a.at["sampleA", "contigs"] == 86387 and a.at["sampleA", "contigs_kept"] == 16094
assert a.at["sampleA", "pct_contigs_kept"] == round(16094 / 86387 * 100, 2)
assert a.at["sampleA", "assembly_length_mb"] == 129.856 and a.at["sampleA", "assembly_mode"] == "assembly"
# reads on contigs and genes from the featureCounts summary
fc = pd.read_csv(f"{stg}/assembly/07_functional_annotation/gene_abundance/sampleA/sampleA.featureCounts.txt.summary",
                 sep="\t", index_col=0).iloc[:, 0]
u = t("st_asm_out", "assembly_read_use").set_index("sample")
assert u.at["sampleA", "pct_reads_on_genes"] == round(fc["Assigned"] / fc.sum() * 100, 2)
assert u.at["sampleA", "pct_reads_on_contigs"] == round((1 - fc["Unassigned_Unmapped"] / fc.sum()) * 100, 2)
qc = t("st_asm_out", "samples_qc").set_index("sample")
assert qc.columns[-3:].tolist() == ["pct_reads_on_genes", "func_contig_genes_annotated_pct",
                                    "func_contig_annotated_by_abundance_pct"], qc.columns[-3:]
assert float(qc.at["sampleA", "pct_reads_on_genes"]) == u.at["sampleA", "pct_reads_on_genes"]
# MIMAG: bin10 has 23S + 16S + 5S and >= 18 tRNA amino acids in Bakta; bin2 (100 % / 0.44 %) has no Bakta, so at
# most medium; the co-assembly's bin17 likewise; sampleB's assembly has no MAG
m = t("st_merge", "mags").set_index(["batch", "assembly", "bin"])
assert m.loc[("asm", "sampleA", "binette_bin10"), "quality_tier"] == "high"
assert m.loc[("asm", "sampleA", "binette_bin10"), "tier_basis"] == "checkm2_rrna_trna"
assert m.loc[("asm", "sampleA", "binette_bin2"), ["quality_tier", "tier_basis"]].tolist() == ["medium", "checkm2_only"]
assert m.loc[("co", "coassembly", "binette_bin17"), "tier_basis"] == "checkm2_only"
assert m.loc[("co", "coassembly", "binette_bin1"), "quality_tier"] == "high"
assert m.loc[("asm", "sampleA", "binette_bin10"), "genus"] == "Senimuribacter"
assert int(m.loc[("asm", "sampleA", "binette_bin10"), "cds"]) == 2291
assert len(m) == 4 and "sampleB" not in m.index.get_level_values("assembly")
# MAG depth: per-sample assembly = length-weighted contig coverage of the bin; co-assembly = the bin depth file
# (sampleC's has the padded fields BEDTOOLS_BIN_DEPTH wrote before)
ct = pd.read_csv(f"{stg}/assembly/05_arg_prediction/contig_level/summary/Contig_tax_and_arg_prediction.tsv",
                 sep="\t", keep_default_na=False)
ct = ct[(ct["Sample"] == "sampleA") & (ct["bin"] == "binette_bin10")]
want = round((ct["coverage"] * ct["contig_length"]).sum() / ct["contig_length"].sum(), 4)
ab = t("st_merge", "mag_abundance").set_index(["batch", "bin", "sample"])
assert ab.at[("asm", "binette_bin10", "sampleA"), "mean_depth_x"] == want, (ab, want)
assert ab.at[("co", "binette_bin1", "sampleC"), "mean_depth_x"] == 0.138756
assert ab.at[("co", "binette_bin1", "sampleD"), "mean_depth_x"] == 0.0519836
ags = pd.read_csv(f"{stg}/coassembly/07_functional_annotation/microbecensus/sampleC/sampleC.ags.tsv", sep="\t")
size = m.loc[("co", "coassembly", "binette_bin1"), "genome_size_mb"] * 1e6
assert ab.at[("co", "binette_bin1", "sampleC"), "pct_sample_bases"] == round(
    0.138756 * size / ags.at[0, "total_bases"] * 100, 4)
# contig ARGs: one row per ARG hit, the bin's tier, the host at the lowest rank Metabuli reached
cl = t("st_asm_out", "arg_contig_long")
src = pd.read_csv(f"{stg}/assembly/05_arg_prediction/contig_level/summary/Contig_tax_and_arg_prediction.tsv",
                  sep="\t", keep_default_na=False)
assert len(cl) == (src["ARG_gene"] != "").sum()
assert set(cl.loc[cl["bin"] == "binette_bin10", "mag_tier"]) == {"high"}
assert set(cl.loc[cl["bin"] == "", "mag_tier"]) <= {""}
assert ((cl["host_rank"] == "unclassified") == (cl["host_taxon"] == "")).all()
h = t("st_asm_out", "arg_contig_hosts")
assert h["arg_genes"].sum() == len(cl)
# a species call counts at its genus
assert "species" not in set(h["host_rank"]) or (cl.loc[cl["host_rank"] == "species", "genus"] == "").any()
# ARG gene abundance: per-sample assemblies repeat gene ids; sampleB's row with a sampleA gene id is not joined
ga = t("st_merge", "arg_contig_abundance")
assert 999 not in set(ga["count"]), ga[ga["count"] == 999]
assert set(ga.loc[ga["assembly"] == "sampleA", "sample"]) == {"sampleA"}
co = ga[ga["assembly"] == "coassembly"]
assert set(co["sample"]) == {"sampleC", "sampleD"} and len(co) == 2 * co["gene_id"].nunique()
# bin ARG clusters: names with underscores, a cluster shared by sampleA's two MAGs
ml = t("st_asm_out", "arg_mag_long")
r = ml[ml["gene_id"] == "k141_8786_20"].iloc[0]
assert (r["arg_gene"], r["arg_class"], r["bin"]) == ("TRUNCATED_PUTATIVE_RESPONSE_REGULATOR_ARLR", "unclassified",
                                                     "binette_bin10"), r
sh = t("st_asm_out", "arg_cluster_sharing")
assert len(sh) == 1 and sh.at[0, "n_mags"] == 2 and sh.at[0, "mags"] == "sampleA:binette_bin10;sampleA:binette_bin2"
# contig functions: eggNOG CAZy names split from the accession, dbCAN apart, COG named; CPGE (MicrobeCensus ran)
cz = pd.read_csv("st_asm_out/tables/function_contig_cazy.tsv", sep="\t", index_col=0, keep_default_na=False)
assert not any("|" in a for a in cz.index) and (cz["description"] != "").all(), cz.head()
assert "function_contig_cazy_dbcan" in open("st_asm_out/README.md").read()
cog = pd.read_csv("st_asm_out/tables/function_contig_cog.tsv", sep="\t", index_col=0)
assert cog.at["S", "description"] == "Function unknown"
fa = pd.read_csv(f"{stg}/assembly/07_functional_annotation/summary/function_abundance.tsv", sep="\t")
row = fa[(fa["ontology"] == "cog") & (fa["accession"] == "K") & (fa["sample_id"] == "sampleA")].iloc[0]
assert abs(cog.at["K", "sampleA"] - row["abundance_cpge"]) < 1e-9
fl = t("st_asm_out", "function_contig_long")
assert len(fl) == len(fa) and "abundance_native" not in fl.columns
assert set(fl.loc[fl["backend"] == "run_dbcan", "ontology"]) == {"cazy_dbcan"}
PYEOF
check "ARG RPKG = RPKM x clean pairs / 1e6 / genome equivalents (NA_no_ags without MicrobeCensus)" py - <<'PYEOF'
import pandas as pd
rpkm = pd.read_csv("st_asm_out/tables/arg_gene_rpkm.tsv", sep="\t", index_col=0, keep_default_na=False)
rpkg = pd.read_csv("st_asm_out/tables/arg_gene_rpkg.tsv", sep="\t", index_col=0, keep_default_na=False)
qc = pd.read_csv("st_asm_out/tables/samples_qc.tsv", sep="\t", index_col=0)
f = qc.at["sampleA", "clean_pairs"] / 1e6 / qc.at["sampleA", "mc_genome_equivalents"]
num = pd.to_numeric(rpkm["sampleA"], errors="coerce").dropna()
assert len(num) and ((pd.to_numeric(rpkg.loc[num.index, "sampleA"]) - num * f).abs() < 1e-9).all()
# run1: sampleB has no MicrobeCensus estimate (and is low depth: its codes stay)
r1 = pd.read_csv("out/tables/arg_gene_rpkg.tsv", sep="\t", index_col=0, keep_default_na=False)
assert set(r1["sampleB"]) <= {"NA_no_ags", "NA_low_depth"}, set(r1["sampleB"])
PYEOF
check "every process of summary_catalog.py exists in the pipeline" py - "${REPO_DIR}" <<'PYEOF'
import re, sys
from pathlib import Path
sys.dont_write_bytecode = True
repo = Path(sys.argv[1])
sys.path.insert(0, str(repo / "bin"))
from summary_catalog import STAGES
text = "\n".join(p.read_text() for p in [repo / "main.nf", *repo.glob("subworkflows/local/*.nf"),
                                         *repo.glob("modules/**/main.nf")])
known = set(re.findall(r"^\s*process\s+(\w+)", text, re.M)) | set(re.findall(r"\bas\s+(\w+)", text))
missing = [p for st in STAGES for p in st["processes"] if p not in known]
assert not missing, missing
PYEOF

echo
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
