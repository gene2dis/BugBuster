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
#   - CARD snps.txt check of variant-model mutations, gene names and tier rules (cfr genes are tier 1)
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

check "CARD snps.txt check, gene names and tier rules" py - "${REPO_DIR}/bin" <<'PYEOF'
import sys, tempfile, os
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
import summary_args as rs
from summary_common import strip_gtdb
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
PYEOF

echo
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
