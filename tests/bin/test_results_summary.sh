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
}

echo "=== results_summary.py test suite ==="
[ "${RS_TEST_HOST:-0}" = "1" ] && echo "python: host" || echo "image: ${SUMMARY_IMG}"
cd "${TMP_DIR}"
make_results run1/results sampleA sampleB

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
PYEOF

echo
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
