#!/bin/bash
#
# Regression tests for audit item #13 (report steps must not mask failures):
#   - bin/report_unify.py exits non-zero when no fastp reports are present
#     (previously exit 0 with no output, hidden by the module's path("*")),
#   - bin/bin_quality_report.py and bin/bin_tax_report.py exit non-zero when
#     every input is header-only (previously a "successful" empty report),
#   - bin/taxonomy_report.py exits non-zero on a missing reads report
#     (previously it silently rebuilt a QC-column-less Reads_report.csv),
#   - and the happy paths still produce their real outputs.
#
# The scripts run inside the same pinned images the report processes use.
#
# Requirements: docker (present on GitHub ubuntu-latest runners).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
FIXTURES="${REPO_DIR}/tests/data/taxonomy"

# Use the exact images pinned in the modules so the test cannot drift from them
MULLED_IMG=$(grep -o "'quay.io/biocontainers/mulled-v2-[^']*'" "${REPO_DIR}/modules/local/reads_report/main.nf" \
    | tr -d "'" | head -1)
REPORT_IMG=$(grep -o "'[^']*community.wave.seqera.io/library/python_pandas[^']*'" "${REPO_DIR}/modules/local/taxonomy_report/main.nf" \
    | tr -d "'" | grep -v '^oras://' | head -1)

if [ -z "${MULLED_IMG}" ] || [ -z "${REPORT_IMG}" ]; then
    echo "ERROR: could not extract pinned container images from the modules"
    exit 1
fi

echo "=== Report script test suite ==="
echo "Mulled image: ${MULLED_IMG}"
echo "Report image: ${REPORT_IMG}"
echo ""

TMP_DIR=$(mktemp -d)
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

WORK="${TMP_DIR}/work"
mkdir -p "${WORK}"

#
# Helpers
#
PASS=0
FAIL=0
LAST_WORKDIR=""
NEXT_WORKDIR=""
RUN_IMG="${MULLED_IMG}"

run_script() {
    # run_script <workdir> <script-name> [args...] — inside ${RUN_IMG}
    local workdir="$1" script="$2"
    shift 2
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e MPLCONFIGDIR=/tmp \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" \
        -v "${workdir}:/repwork" -w /repwork \
        "${RUN_IMG}" python3 "/pipeline_bin/${script}" "$@"
}

expect_pass() {
    # honors NEXT_WORKDIR when the caller pre-seeded a directory with fixtures
    local desc="$1"; shift
    local workdir="${NEXT_WORKDIR:-${WORK}/$(echo "${desc}" | tr ' /:' '___')}"
    NEXT_WORKDIR=""
    mkdir -p "${workdir}"
    LAST_WORKDIR="${workdir}"
    if run_script "${workdir}" "$@" > "${workdir}.log" 2>&1; then
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${desc} — expected success, got failure:"
        tail -5 "${workdir}.log" | sed 's/^/    /'
        FAIL=$((FAIL + 1))
    fi
}

expect_fail() {
    local desc="$1"; shift
    local workdir="${NEXT_WORKDIR:-${WORK}/$(echo "${desc}" | tr ' /:' '___')}"
    NEXT_WORKDIR=""
    mkdir -p "${workdir}"
    LAST_WORKDIR="${workdir}"
    if run_script "${workdir}" "$@" > "${workdir}.log" 2>&1; then
        echo "✗ ${desc} — expected nonzero exit, script succeeded"
        FAIL=$((FAIL + 1))
    else
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    fi
}

check_grep() {
    # check_grep <desc> <path-relative-to-last-workdir> <pattern>
    if grep -q "$3" "${LAST_WORKDIR}/$2" 2>/dev/null; then
        echo "✓ $1"
        PASS=$((PASS + 1))
    else
        echo "✗ $1 — pattern '$3' not found in $2"
        FAIL=$((FAIL + 1))
    fi
}

check_file() {
    # check_file <desc> <path-relative-to-last-workdir>
    if [ -e "${LAST_WORKDIR}/$2" ]; then
        echo "✓ $1"
        PASS=$((PASS + 1))
    else
        echo "✗ $1 — missing: $2"
        FAIL=$((FAIL + 1))
    fi
}

#
# report_unify.py (fixture formats match modules/local/qfilter,
# bowtie2_decontaminate, and count_reads report outputs)
#
echo "--- report_unify.py ---"
RU_WORK="${WORK}/report_unify_filtered"
mkdir -p "${RU_WORK}"
printf 'Id\tRaw reads\tFastp\ns1\t1000\t900\n' > "${RU_WORK}/s1_fastp_report.tsv"
printf 'Id\tClean reads (contaminants removed)\ns1\t850\n' > "${RU_WORK}/s1_contaminants_decontamination_report.tsv"
NEXT_WORKDIR="${RU_WORK}"
expect_pass "report_unify: filtered mode with fastp + bowtie reports" report_unify.py contaminants
check_grep "report_unify: merged report header" "Reads_report.csv" \
    "Id,Raw reads,Fastp,Clean reads (contaminants removed)"
check_grep "report_unify: merged report row" "Reads_report.csv" "^s1,1000,900,850$"
check_file "report_unify: boxplot created" "Box_plot_reads.png"

RU_NONE="${WORK}/report_unify_none"
mkdir -p "${RU_NONE}"
printf 'Id\tRaw reads\ns1\t1000\n' > "${RU_NONE}/s1_fastp_report.tsv"
NEXT_WORKDIR="${RU_NONE}"
expect_pass "report_unify: none mode with count report" report_unify.py none
check_grep "report_unify: none-mode report row" "Reads_report.csv" "^s1,1000$"

# No fastp reports staged: both modes must exit non-zero (the masking bug)
RU_EMPTY="${WORK}/report_unify_empty_filtered"
mkdir -p "${RU_EMPTY}"
printf 'Id\tClean reads (contaminants removed)\ns1\t850\n' > "${RU_EMPTY}/s1_contaminants_decontamination_report.tsv"
NEXT_WORKDIR="${RU_EMPTY}"
expect_fail "report_unify: filtered mode without fastp reports fails" report_unify.py contaminants
expect_fail "report_unify: none mode without fastp reports fails" report_unify.py none

#
# taxonomy_report.py — missing required reads report must fail (previously it
# silently rebuilt a QC-column-less Reads_report.csv). Runs in the report image.
#
echo "--- taxonomy_report.py ---"
RUN_IMG="${REPORT_IMG}"
TR_WORK="${WORK}/taxonomy_report_missing_reads"
mkdir -p "${TR_WORK}"
cp "${FIXTURES}/sample1.kraken2.report.txt" "${TR_WORK}/"
NEXT_WORKDIR="${TR_WORK}"
expect_fail "taxonomy_report: missing reads report fails" taxonomy_report.py --profiler kraken2 \
    --reports sample1.kraken2.report.txt \
    --reads-report nonexistent_reads_report.csv --db-name standard-8 --output-dir .
check_grep "taxonomy_report: failure names the missing input" "../taxonomy_report_missing_reads.log" \
    "Required reads report"
RUN_IMG="${MULLED_IMG}"

#
# bin_quality_report.py (CheckM2 batch report format)
#
echo "--- bin_quality_report.py ---"
BQ_WORK="${WORK}/bin_quality_happy"
mkdir -p "${BQ_WORK}"
printf 'Name\tCompleteness\tContamination\nbin.1\t95.0\t2.0\nbin.2\t60.0\t4.0\n' \
    > "${BQ_WORK}/sampleA_semibin_quality_report.tsv"
printf 'Name\tCompleteness\tContamination\nbin.1\t96.0\t1.0\n' \
    > "${BQ_WORK}/sampleA_metawrap_quality_report.tsv"
NEXT_WORKDIR="${BQ_WORK}"
expect_pass "bin_quality: real CheckM2 reports summarized" bin_quality_report.py
check_grep "bin_quality: summary counts High semibin MAG" "Mag_quality_summary.csv" "^Semibin,High,1$"
check_grep "bin_quality: all-MAGs table row" "All_Mag_quality_table.csv" "^sampleA,Semibin,bin.1,High,95.0,2.0$"
check_grep "bin_quality: refined table has metawrap MAG" "Refined_Mag_quality_table.csv" "^sampleA,Metawrap,bin.1,High,96.0,1.0$"
check_file "bin_quality: quality plot created" "Total_bins_quality_plot.png"

BQ_EMPTY="${WORK}/bin_quality_empty"
mkdir -p "${BQ_EMPTY}"
printf 'Name\tCompleteness\tContamination\n' > "${BQ_EMPTY}/sampleA_semibin_quality_report.tsv"
NEXT_WORKDIR="${BQ_EMPTY}"
expect_fail "bin_quality: all-header-only input fails" bin_quality_report.py
check_grep "bin_quality: failure message names the empty input" "../bin_quality_empty.log" "no bins were assessed"

#
# bin_tax_report.py (GTDB-Tk batch summary format)
#
echo "--- bin_tax_report.py ---"
BT_WORK="${WORK}/bin_tax_happy"
mkdir -p "${BT_WORK}"
printf 'user_genome\tclassification\tfastani_reference\nbin.1\td__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis\tGCF_000009045.1\n' \
    > "${BT_WORK}/sampleA_gtdbtk_bac120.tsv"
NEXT_WORKDIR="${BT_WORK}"
expect_pass "bin_tax: real GTDB-Tk summary parsed" bin_tax_report.py
check_grep "bin_tax: MAG row with parsed ranks" "MAGs_tax_summary.csv" \
    "^sampleA,bin.1,GCF_000009045.1,Bacteria,Bacillota,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus subtilis$"
check_file "bin_tax: taxonomy plot created" "MAGs_tax_plot.png"

BT_EMPTY="${WORK}/bin_tax_empty"
mkdir -p "${BT_EMPTY}"
printf 'user_genome\tclassification\tfastani_reference\n' > "${BT_EMPTY}/sampleA_gtdbtk_bac120.tsv"
NEXT_WORKDIR="${BT_EMPTY}"
expect_fail "bin_tax: all-header-only input fails" bin_tax_report.py
check_grep "bin_tax: failure message names the empty input" "../bin_tax_empty.log" "no bins were classified"

#
# bin_summary.py — audit item #11: ALL staged quality and taxonomy files must be
# read (previously glob(...)[0] kept one arbitrary file of each, silently
# dropping other binners' bins and whole domains such as archaea from ar53).
#
echo "--- bin_summary.py ---"
BS_WORK="${WORK}/bin_summary_happy"
mkdir -p "${BS_WORK}"
# Depth files as produced by BEDTOOLS (one per read sample in co-assembly mode)
printf 'sample\tbin_id\tTotal_contigs\taverage_cov\ns1\tbin.1\t10\t12.5\ns1\tbin.2\t8\t3.0\ns1\tbin.3\t5\t7.5\n' \
    > "${BS_WORK}/s1_bin_depth.tsv"
printf 'sample\tbin_id\tTotal_contigs\taverage_cov\ns2\tbin.1\t10\t2.5\ns2\tbin.2\t8\t9.0\ns2\tbin.3\t5\t1.5\n' \
    > "${BS_WORK}/s2_bin_depth.tsv"
# TWO quality reports: bin.3 lives only in the second one
printf 'Name\tCompleteness\tContamination\tContig_N50\tGenome_Size\tTotal_Coding_Sequences\nbin.1\t95.0\t2.0\t50000\t3000000\t2900\nbin.2\t90.0\t1.0\t40000\t2000000\t2100\n' \
    > "${BS_WORK}/coassembly_semibin_quality_report.tsv"
printf 'Name\tCompleteness\tContamination\tContig_N50\tGenome_Size\tTotal_Coding_Sequences\nbin.3\t85.0\t3.0\t30000\t2500000\t2400\n' \
    > "${BS_WORK}/coassembly_metabat_quality_report.tsv"
# bac120 + ar53 summaries: the archaeal bin lives only in ar53; plus a zero-byte
# ar53 file (what GTDB_TK_BATCH stubs/no-archaea samples produce) that must be skipped
printf 'user_genome\tclassification\nbin.1\td__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis\nbin.3\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\n' \
    > "${BS_WORK}/coassembly_gtdbtk_bac120.tsv"
printf 'user_genome\tclassification\nbin.2\td__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter smithii\n' \
    > "${BS_WORK}/coassembly_gtdbtk_ar53.tsv"
: > "${BS_WORK}/empty_gtdbtk_ar53.tsv"
NEXT_WORKDIR="${BS_WORK}"
expect_pass "bin_summary: all quality and taxonomy files combined" bin_summary.py
check_grep "bin_summary: bin from second quality report present" "Bin_summary.csv" "^bin.3,85.0,3.0"
check_grep "bin_summary: archaeal bin from ar53 present" "Bin_summary.csv" "^bin.2,.*Archaea,Methanobacteriota"
check_grep "bin_summary: per-sample coverage columns pivoted" "Bin_summary.csv" "s1_avgcov,s2_avgcov"
BS_ROWS=$(tail -n +2 "${BS_WORK}/Bin_summary.csv" 2>/dev/null | wc -l)
if [ "${BS_ROWS}" -eq 3 ]; then
    echo "✓ bin_summary: summary has all 3 bins"
    PASS=$((PASS + 1))
else
    echo "✗ bin_summary: expected 3 rows, got ${BS_ROWS}"
    FAIL=$((FAIL + 1))
fi

# All quality reports header-only: must exit non-zero, not write an empty summary
BS_EMPTY="${WORK}/bin_summary_empty"
mkdir -p "${BS_EMPTY}"
cp "${BS_WORK}/s1_bin_depth.tsv" "${BS_WORK}/coassembly_gtdbtk_bac120.tsv" "${BS_EMPTY}/"
printf 'Name\tCompleteness\tContamination\tContig_N50\tGenome_Size\tTotal_Coding_Sequences\n' \
    > "${BS_EMPTY}/coassembly_semibin_quality_report.tsv"
NEXT_WORKDIR="${BS_EMPTY}"
expect_fail "bin_summary: all-header-only quality input fails" bin_summary.py

#
# arg_norm_report.py — audit item #22: a sample where ARGs-OAP finds nothing
# (nRead/nCell = 0) must yield empty normalization cells, not inf, in the
# published CSVs.
#
echo "--- arg_norm_report.py ---"
AN_WORK="${WORK}/arg_norm_zero_cells"
mkdir -p "${AN_WORK}"
KARGA_GENE='MEG_1|Drugs|Aminoglycosides|Aminoglycoside_resistance|APH3'
KARGVA_GENE='ID0|K12A|x|Escherichia coli resistance to betalactams|x|BETALACTAM|Betalactams|multi|blaTEM|x|x|x|x|blaTEM_1|NA|x|E_coli_resistant|end'
for s in s1 s2; do
    printf 'GeneIdx,PercentGeneCovered,AverageKMerDepth\n%s,95%%,10.0\n' "${KARGA_GENE}" \
        > "${AN_WORK}/${s}_all_reads_KARGA_mappedGenes.csv"
    printf 'GeneIdx,KmerSNPHits,PercentGeneCovered,AverageKMerDepth\n%s,3/10,85%%,10.0\n' "${KARGVA_GENE}" \
        > "${AN_WORK}/${s}_all_reads_KARGVA_mappedGenes.csv"
done
mkdir -p "${AN_WORK}/s1_args_oap_s1_out" "${AN_WORK}/s2_args_oap_s1_out"
printf 'sample\tnRead\tn16S\tnCell\ns1\t1000000\t50\t100\n' > "${AN_WORK}/s1_args_oap_s1_out/metadata.txt"
# s2: ARGs-OAP found no reads/16S — the division-by-zero edge case
printf 'sample\tnRead\tn16S\tnCell\ns2\t0\t0\t0\n' > "${AN_WORK}/s2_args_oap_s1_out/metadata.txt"
NEXT_WORKDIR="${AN_WORK}"
expect_pass "arg_norm: zero-count sample handled" arg_norm_report.py
check_grep "arg_norm: KARGA normalization computed for s1" "KARGA_norm.csv" "^s1,1000000,50,100,.*,10.0,0.1$"
check_file "arg_norm: KARGVA report created" "KARGVA_norm.csv"
for f in KARGA_norm.csv KARGVA_norm.csv; do
    if grep -qi "inf" "${AN_WORK}/${f}" 2>/dev/null; then
        echo "✗ arg_norm: ${f} contains inf"
        FAIL=$((FAIL + 1))
    else
        echo "✓ arg_norm: ${f} has no inf"
        PASS=$((PASS + 1))
    fi
done

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
