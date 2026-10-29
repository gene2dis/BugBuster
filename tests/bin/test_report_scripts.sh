#!/bin/bash
#
# Regression tests for audit item #13 (report steps must not mask failures):
#   - bin/report_unify.py exits non-zero when no fastp reports are present
#     (previously exit 0 with no output, hidden by the module's path("*")),
#   - bin/bin_quality_report.py, bin/bin_tax_report.py and bin/bin_summary.py
#     warn and write empty outputs when every input is header-only (SemiBin
#     and Binette finish with empty bin sets, so a run without bins
#     finishes), but still fail when no input file is staged,
#   - bin/taxonomy_report.py exits non-zero on a missing reads report
#     (previously it silently rebuilt a QC-column-less Reads_report.csv),
#   - and the happy paths still produce their real outputs,
#   - bin/arg_contig_level_report.py writes a real TSV, and bin/arg_blobplot.py
#     plots only the ARG-carrying contigs (contigs without an ARG were "Other"),
#   - bin/blobplot.py and bin/arg_blobplot.py pickle their figure without
#     error (design doc T9 finding: lambda tick formatters broke pickle.dump).
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
# REPORT_TEST_IMG=<image>: run the scripts in that image instead of the pinned one (to try an image before pinning it)
MULLED_IMG=${REPORT_TEST_IMG:-${MULLED_IMG}}
REPORT_IMG=${REPORT_TEST_IMG:-${REPORT_IMG}}

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
    > "${BQ_WORK}/sampleA_binette_quality_report.tsv"
NEXT_WORKDIR="${BQ_WORK}"
expect_pass "bin_quality: real CheckM2 reports summarized" bin_quality_report.py
check_grep "bin_quality: summary counts High semibin MAG" "Mag_quality_summary.csv" "^Semibin,High,1$"
check_grep "bin_quality: all-MAGs table row" "All_Mag_quality_table.csv" "^sampleA,Semibin,bin.1,High,95.0,2.0$"
check_grep "bin_quality: refined table has Binette MAG" "Refined_Mag_quality_table.csv" "^sampleA,Binette,bin.1,High,96.0,1.0$"
check_file "bin_quality: quality plot created" "Total_bins_quality_plot.png"

BQ_EMPTY="${WORK}/bin_quality_empty"
mkdir -p "${BQ_EMPTY}"
printf 'Name\tCompleteness\tContamination\n' > "${BQ_EMPTY}/sampleA_semibin_quality_report.tsv"
NEXT_WORKDIR="${BQ_EMPTY}"
# No sample with a bin: warn and write empty outputs
expect_pass "bin_quality: all-header-only input warns and succeeds" bin_quality_report.py
check_grep "bin_quality: warning names the empty input" "../bin_quality_empty.log" "WARNING: all CheckM2 quality reports are empty"
check_grep "bin_quality: header-only summary written" "Mag_quality_summary.csv" "^Binner,Quality classify,Total MAGs$"
check_grep "bin_quality: header-only all-MAGs table written" "All_Mag_quality_table.csv" \
    "^Sample,Binner,Bin name,Quality classify,Completeness,Contamination$"
check_grep "bin_quality: header-only refined table written" "Refined_Mag_quality_table.csv" \
    "^Sample,Binner,Bin name,Quality classify,Completeness,Contamination$"
check_file "bin_quality: placeholder plot written" "Total_bins_quality_plot.png"

BQ_NONE="${WORK}/bin_quality_none"
mkdir -p "${BQ_NONE}"
NEXT_WORKDIR="${BQ_NONE}"
expect_fail "bin_quality: no CheckM2 file still fails" bin_quality_report.py

#
# bin_tax_report.py (GTDB-Tk batch summary format)
#
echo "--- bin_tax_report.py ---"
# GTDB-Tk >= 2.7 (the pinned 2.8.0) names the reference column
# closest_genome_reference; the script's fallback also accepts the legacy
# fastani_reference name (<= 2.6-era summaries). Both branches are pinned.
BT_WORK="${WORK}/bin_tax_happy"
mkdir -p "${BT_WORK}"
printf 'user_genome\tclassification\tclosest_genome_reference\nbin.1\td__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis\tGCF_000009045.1\n' \
    > "${BT_WORK}/sampleA_gtdbtk_bac120.tsv"
NEXT_WORKDIR="${BT_WORK}"
expect_pass "bin_tax: real GTDB-Tk 2.7 summary parsed" bin_tax_report.py
check_grep "bin_tax: MAG row with parsed ranks" "MAGs_tax_summary.csv" \
    "^sampleA,bin.1,GCF_000009045.1,Bacteria,Bacillota,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus subtilis$"
check_file "bin_tax: taxonomy plot created" "MAGs_tax_plot.png"

BT_LEGACY="${WORK}/bin_tax_legacy"
mkdir -p "${BT_LEGACY}"
printf 'user_genome\tclassification\tfastani_reference\nbin.1\td__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis\tGCF_000009045.1\n' \
    > "${BT_LEGACY}/sampleA_gtdbtk_bac120.tsv"
NEXT_WORKDIR="${BT_LEGACY}"
expect_pass "bin_tax: legacy fastani_reference column still parsed" bin_tax_report.py
check_grep "bin_tax: legacy summary MAG row" "MAGs_tax_summary.csv" \
    "^sampleA,bin.1,GCF_000009045.1,Bacteria,"

BT_EMPTY="${WORK}/bin_tax_empty"
mkdir -p "${BT_EMPTY}"
printf 'user_genome\tclassification\tclosest_genome_reference\n' > "${BT_EMPTY}/sampleA_gtdbtk_bac120.tsv"
NEXT_WORKDIR="${BT_EMPTY}"
# No sample with a classified bin: warn and write empty outputs (SemiBin and
# Binette finish with empty bin sets, so a binless run must finish too)
expect_pass "bin_tax: all-header-only input warns and succeeds" bin_tax_report.py
check_grep "bin_tax: warning names the empty input" "../bin_tax_empty.log" "WARNING: all GTDB-Tk summaries are empty"
check_grep "bin_tax: header-only summary written" "MAGs_tax_summary.csv" \
    "^sample,bin_name,closest_genome_reference,Domain,Phylum,Class,Order,Family,Genus,Species$"
check_file "bin_tax: placeholder plot written" "MAGs_tax_plot.png"

BT_NONE="${WORK}/bin_tax_none"
mkdir -p "${BT_NONE}"
NEXT_WORKDIR="${BT_NONE}"
expect_fail "bin_tax: no GTDB-Tk file still fails" bin_tax_report.py

#
# bin_summary.py — audit item #11: ALL staged quality and taxonomy files must be
# read (previously glob(...)[0] kept one arbitrary file of each, silently
# dropping other binners' bins and whole domains such as archaea from ar53).
#
echo "--- bin_summary.py ---"
BS_WORK="${WORK}/bin_summary_happy"
mkdir -p "${BS_WORK}"
# Depth files as produced by BEDTOOLS_BIN_DEPTH (one per read sample in co-assembly mode)
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

# Quality reports header-only but classified bins present: CheckM2 lost
# bins, so it must exit non-zero, not write an empty summary
BS_EMPTY="${WORK}/bin_summary_empty"
mkdir -p "${BS_EMPTY}"
cp "${BS_WORK}/s1_bin_depth.tsv" "${BS_WORK}/coassembly_gtdbtk_bac120.tsv" "${BS_EMPTY}/"
printf 'Name\tCompleteness\tContamination\tContig_N50\tGenome_Size\tTotal_Coding_Sequences\n' \
    > "${BS_EMPTY}/coassembly_semibin_quality_report.tsv"
NEXT_WORKDIR="${BS_EMPTY}"
expect_fail "bin_summary: header-only quality with classified bins fails" bin_summary.py

# No refined bin at all: header-only quality and taxonomy, no depth files
# (the bin depth step is skipped). Warn and write a header-only summary
BS_NOBINS="${WORK}/bin_summary_nobins"
mkdir -p "${BS_NOBINS}"
printf 'Name\tCompleteness\tContamination\tContig_N50\tGenome_Size\tTotal_Coding_Sequences\n' \
    > "${BS_NOBINS}/coassembly_binette_quality_report.tsv"
printf 'user_genome\tclassification\n' > "${BS_NOBINS}/coassembly_gtdbtk_bac120.tsv"
: > "${BS_NOBINS}/coassembly_gtdbtk_ar53.tsv"
NEXT_WORKDIR="${BS_NOBINS}"
expect_pass "bin_summary: no refined bins warns and succeeds" bin_summary.py
check_grep "bin_summary: warning names the empty input" "../bin_summary_nobins.log" "WARNING: all quality reports and GTDB-Tk summaries are empty"
check_grep "bin_summary: header-only summary written" "Bin_summary.csv" \
    "^Name,Completeness,Contamination,Contig_N50,Genome_Size,Total_Coding_Sequences,Total_contigs,Domain,Phylum,Class,Order,Family,Genus,Species$"

BS_NONE="${WORK}/bin_summary_none"
mkdir -p "${BS_NONE}"
printf 'user_genome\tclassification\n' > "${BS_NONE}/coassembly_gtdbtk_bac120.tsv"
NEXT_WORKDIR="${BS_NONE}"
expect_fail "bin_summary: no quality file still fails" bin_summary.py

#
# kmarvel_report.py — real K-MARVEL tables (tests/data/kmarvel: sampleA with
# two ARGs, sampleB header-only); the zero-ARG sample must stay as an
# all-zero column. Runs in the mulled image KMARVEL_REPORT pins (MULLED_IMG).
#
echo "--- kmarvel_report.py ---"
KM_WORK="${WORK}/kmarvel_report"
mkdir -p "${KM_WORK}"
cp "${REPO_DIR}"/tests/data/kmarvel/sample{A,B}_kmarvel.tsv "${KM_WORK}/"
NEXT_WORKDIR="${KM_WORK}"
expect_pass "kmarvel_report: ARG and zero-ARG samples merged" kmarvel_report.py \
    sampleA_kmarvel.tsv sampleB_kmarvel.tsv
check_grep "kmarvel_report: matrix keeps the zero-ARG sample column" \
    "KMARVEL_supporting_reads_matrix.tsv" "^aro_id	aro_full_name	sampleA	sampleB$"
check_grep "kmarvel_report: supporting reads per ARG and sample" \
    "KMARVEL_supporting_reads_matrix.tsv" "^ARO:3000980	TEM-117	60	0$"
check_grep "kmarvel_report: read pairs summed per AMR gene family" \
    "KMARVEL_gene_family_matrix.tsv" "^TEM beta-lactamase	60	0$"
check_grep "kmarvel_report: multi-class ARG counted once per drug class" \
    "KMARVEL_drug_class_matrix.tsv" "^cephalosporin	2	0$"
check_file "kmarvel_report: long table created" "KMARVEL_summary_long.tsv"
NEXT_WORKDIR="${KM_WORK}"
expect_pass "kmarvel_report: only zero-ARG samples still reports" kmarvel_report.py sampleB_kmarvel.tsv
check_grep "kmarvel_report: empty matrix keeps the sample header" \
    "KMARVEL_supporting_reads_matrix.tsv" "^aro_id	aro_full_name	sampleB$"
check_grep "kmarvel_report: empty gene-family matrix keeps the sample header" \
    "KMARVEL_gene_family_matrix.tsv" "^gene_family	sampleB$"
expect_fail "kmarvel_report: no input tables fails" kmarvel_report.py
KM_BAD="${WORK}/kmarvel_report_bad"
mkdir -p "${KM_BAD}"
cut -f1-5 "${REPO_DIR}/tests/data/kmarvel/sampleA_kmarvel.tsv" > "${KM_BAD}/s1_kmarvel.tsv"
NEXT_WORKDIR="${KM_BAD}"
expect_fail "kmarvel_report: table without the K-MARVEL columns fails" kmarvel_report.py s1_kmarvel.tsv

#
# args_oap_report.py — ARGs-OAP 4.0.0 stage two rows (tests/data/args_oap:
# sampleA with 12 SARG hits, two from its singletons; sampleB metadata.txt
# only). A sample is its pairs plus singletons, normalised over both: tetA(48)
# 4 reads, 0.548105 copies / 16 cells. Runs in the mulled image ARGS_OAP_REPORT
# pins (MULLED_IMG).
#
echo "--- args_oap_report.py ---"
AO_WORK="${WORK}/args_oap_report"
mkdir -p "${AO_WORK}"
cp -r "${REPO_DIR}"/tests/data/args_oap/sample{A,B}_args_oap "${AO_WORK}/"
mkdir -p "${AO_WORK}/sampleC_args_oap"
# stage two stops without tables when a sample has no 16S or cells
printf 'sample\tnRead\tn16S\tnCell\nsampleC_tmp_reads\t2000\t0.0\t0.0\n' > "${AO_WORK}/sampleC_args_oap/metadata.txt"
NEXT_WORKDIR="${AO_WORK}"
expect_pass "args_oap_report: samples with, without and with unnormalisable SARG hits" args_oap_report.py \
    sampleA_args_oap sampleB_args_oap sampleC_args_oap
check_grep "args_oap_report: pairs and singletons summed per sample" \
    "args_oap_cells.tsv" "^sampleA	120000	32.0	16.0	9.0	ok$"
check_grep "args_oap_report: no SARG hit" "args_oap_cells.tsv" "^sampleB	4000	1.0	0.5	0.0	no_sarg_hit$"
check_grep "args_oap_report: no 16S or cells" "args_oap_cells.tsv" "^sampleC	2000	0.0	0.0	0.0	not_normalised$"
check_grep "args_oap_report: type level over pairs and singletons" \
    "args_oap_long.tsv" "^sampleA	type	tetracycline			4.0	0.548104956268221[0-9]*	0.03425655976676[0-9]*	0.0171282798833[0-9]*	33.33333"
check_grep "args_oap_report: 1/3-weighted reads of a multi-component reference" \
    "args_oap_long.tsv" "^sampleA	subtype	multidrug	multidrug__adeF		1.0	"
check_grep "args_oap_report: gene level keeps type and subtype" \
    "args_oap_long.tsv" "^sampleA	gene	polymyxin	polymyxin__rosB	ZP_04612557	0.5	"
if grep -q "^sample[BC]	" "${AO_WORK}/args_oap_long.tsv"; then
    echo "✗ args_oap_report: samples without hits have no SARG rows"; FAIL=$((FAIL + 1))
else
    echo "✓ args_oap_report: samples without hits have no SARG rows"; PASS=$((PASS + 1))
fi
# ARGS_OAP's stage_two_status.txt decides, not the summed counts: summed pairs + singletons can be non-zero
# while no single ARGs-OAP row had reads, 16S and cells (stage two skipped)
AO_MARK="${WORK}/args_oap_report_marker"
mkdir -p "${AO_MARK}/sampleF_args_oap"
printf 'sample\tnRead\tn16S\tnCell\nsampleF_tmp_reads\t2000\t0.0\t0.4\nsampleF_tmp_reads_S\t100\t0.1\t0.0\n' \
    > "${AO_MARK}/sampleF_args_oap/metadata.txt"
echo not_normalised > "${AO_MARK}/sampleF_args_oap/stage_two_status.txt"
NEXT_WORKDIR="${AO_MARK}"
expect_pass "args_oap_report: stage two status taken from ARGS_OAP" args_oap_report.py sampleF_args_oap
check_grep "args_oap_report: marker overrides the summed counts" "args_oap_cells.tsv" "^sampleF	2100	0.1	0.4	0.0	not_normalised$"
echo maybe > "${AO_MARK}/sampleF_args_oap/stage_two_status.txt"
NEXT_WORKDIR="${AO_MARK}"
expect_fail "args_oap_report: an unknown stage two status fails" args_oap_report.py sampleF_args_oap
AO_BAD="${WORK}/args_oap_report_bad"
mkdir -p "${AO_BAD}/sampleD_args_oap" "${AO_BAD}/sampleE_out"
NEXT_WORKDIR="${AO_BAD}"
expect_fail "args_oap_report: a directory without metadata.txt fails" args_oap_report.py sampleD_args_oap
NEXT_WORKDIR="${AO_BAD}"
expect_fail "args_oap_report: a directory not named <sample>_args_oap fails" args_oap_report.py sampleE_out

#
# blobplot.py / arg_blobplot.py (design doc T9 finding, 2026-10-09): the
# figure is also pickled, and lambda tick formatters made pickle.dump fail
# after the PNG was written. Fixtures are cut from a real T9 run (BlobTools
# 1.1.1 table, ARG_CONTIG_LEVEL_REPORT output); runs in the image the
# BLOBPLOT/ARG_BLOBPLOT modules pin.
#
echo "--- blobplot.py / arg_blobplot.py ---"
RUN_IMG=$(grep -o "'quay.io/biocontainers/mulled-v2-[^']*'" "${REPO_DIR}/modules/local/blobplot/main.nf" \
    | tr -d "'" | head -1)
check_pickle() {
    # check_pickle <desc> <script-module> <pkl> — the pickled figure loads back
    if docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e MPLCONFIGDIR=/tmp \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" -v "${LAST_WORKDIR}:/repwork" -w /repwork \
        -e PYTHONPATH=/pipeline_bin "${RUN_IMG}" python3 -c \
        "import pickle, sys; sys.modules['__main__'] = __import__('$2'); fig = pickle.load(open('$3', 'rb')); assert len(fig.axes) == 2" \
        > /dev/null 2>&1; then
        echo "✓ $1"
        PASS=$((PASS + 1))
    else
        echo "✗ $1"
        FAIL=$((FAIL + 1))
    fi
}
BP_WORK="${WORK}/blobplot"
mkdir -p "${BP_WORK}"
cp "${REPO_DIR}/tests/data/reports/test_sample2_Blob_table.test_sample2.blobDB.table.txt" "${BP_WORK}/"
NEXT_WORKDIR="${BP_WORK}"
expect_pass "blobplot: real BlobTools table plotted and pickled" blobplot.py
check_file "blobplot: PNG created" "Phylum_blob_plot.png"
check_pickle "blobplot: pickled figure loads back" blobplot "Phylum_final_blobplot.pkl"
ABP_WORK="${WORK}/arg_blobplot"
mkdir -p "${ABP_WORK}"
cp "${REPO_DIR}/tests/data/reports/Contig_tax_and_arg_prediction.tsv" "${ABP_WORK}/"
NEXT_WORKDIR="${ABP_WORK}"
expect_pass "arg_blobplot: contig ARG table plotted and pickled" arg_blobplot.py
check_file "arg_blobplot: PNG created" "ARG_blob_plot.png"
check_pickle "arg_blobplot: pickled figure loads back" arg_blobplot "ARG_final_blobplot.pkl"
# The fixture has 7 ARG-carrying contigs out of 40 in 5 classes; the 33
# contigs without an ARG used to be plotted as "Other"
check_grep "arg_blobplot: only the ARG-carrying contigs are plotted" "../arg_blobplot.log" "Total contigs: 7$"
check_grep "arg_blobplot: five ARG classes and no 'Other'" "../arg_blobplot.log" \
    "Class distribution: {'Aminoglycoside': 14.29, 'Beta-Lactam': 28.57, 'Fosmidomycin': 14.29, 'Multidrug': 28.57, 'Tetracycline': 14.29}$"

#
# arg_contig_level_report.py: the contig table is written tab-separated, as
# its .tsv name says (it was comma-separated). Inputs: the real BlobTools
# table above (length, GC, coverage), a Metabuli classifications file for its
# contigs (lines in Metabuli 1.2.0's --lineage 1 format, trailing tab
# included: species, genus, a virus and an unclassified contig) and a
# synthetic DeepARG mapping file (DeepARG 1.0.4 header, one acrA hit on
# contig k141_0, values taken from a real test_sample1 hit).
#
echo "--- arg_contig_level_report.py ---"
ACR_WORK="${WORK}/arg_contig_level_report"
mkdir -p "${ACR_WORK}"
cp "${REPO_DIR}/tests/data/reports/test_sample2_Blob_table.test_sample2.blobDB.table.txt" \
   "${REPO_DIR}/tests/data/reports/test_sample2_classifications.tsv" \
   "${REPO_DIR}/tests/data/reports/test_sample2_contigs_deep_arg.out.mapping.ARG" "${ACR_WORK}/"
NEXT_WORKDIR="${ACR_WORK}"
expect_pass "arg_contig_level_report: BlobTools, Metabuli and DeepARG tables joined" arg_contig_level_report.py
check_grep "arg_contig_level_report: tab-separated header with Metabuli's seven ranks" "Contig_tax_and_arg_prediction.tsv" \
    "^Sample	contig_name	contig_length	GC	coverage	superkingdom	phylum	class	order	family	genus	species	taxonomy_rank	taxonomy_score	ARG_gene	.*	counts$"
check_grep "arg_contig_level_report: ARG joined to its contig, ranks from the lineage" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_0	1196	0.4005	3.7229	Bacteria	Pseudomonadota	Gammaproteobacteria	Moraxellales	Moraxellaceae	Acinetobacter	Acinetobacter baumannii	species	0.95811	ACRA	multidrug	"
check_grep "arg_contig_level_report: genus-level call leaves species empty" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_820	1015	0.4138	1.9009	Bacteria	.*	Acinetobacter		genus	0.412783							$"
check_grep "arg_contig_level_report: Viruses (acellular root) is the superkingdom" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_231	[^	]*	[^	]*	[^	]*	Viruses	Uroviricota	Caudoviricetes	Crassvirales	Suoliviridae	Cohcovirus		genus	"
check_grep "arg_contig_level_report: unclassified contig kept with empty ranks" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_456	[^	]*	[^	]*	[^	]*								unclassified	0.0	"
ACR_NOTAX="${WORK}/arg_contig_level_report_notax"
mkdir -p "${ACR_NOTAX}"
cp "${REPO_DIR}/tests/data/reports/test_sample2_Blob_table.test_sample2.blobDB.table.txt" \
   "${REPO_DIR}/tests/data/reports/test_sample2_contigs_deep_arg.out.mapping.ARG" "${ACR_NOTAX}/"
NEXT_WORKDIR="${ACR_NOTAX}"
expect_fail "arg_contig_level_report: no Metabuli file fails" arg_contig_level_report.py

#
# bin_contig_labels.py + the report's bin columns (Phase 9b, L5 a): contig ->
# refined bin (from the bins alone); a contig in two bins counts once (first
# bin, warning); a bins dir with only a marker gives a header. The report adds
# the bin's GTDB-Tk string, untouched, from <sample>_gtdbtk_*.tsv; a bin
# GTDB-Tk has no row for keeps an empty classification.
#
echo "--- bin_contig_labels.py ---"
BCL_WORK="${WORK}/bin_contig_labels"
mkdir -p "${BCL_WORK}/test_sample2_binette_bins" "${BCL_WORK}/empty_bins"
printf ">k141_0 flag=1\nACGT\n>k141_820\nACGT\n" > "${BCL_WORK}/test_sample2_binette_bins/binette_bin1.fa"
printf ">k141_231\nACGT\n>k141_0\nACGT\n" > "${BCL_WORK}/test_sample2_binette_bins/binette_bin2.fa"
touch "${BCL_WORK}/empty_bins/SKIPPED.txt"
printf "user_genome\tclassification\tclosest_genome_reference\nbinette_bin1\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__Acinetobacter baumannii\tGCF_009035845.1\n" \
    > "${BCL_WORK}/test_sample2_gtdbtk_bac120.tsv"
printf "user_genome\tclassification\n" > "${BCL_WORK}/test_sample2_gtdbtk_ar53.tsv"
NEXT_WORKDIR="${BCL_WORK}"
expect_pass "bin_contig_labels: contig to bin table" bin_contig_labels.py test_sample2 test_sample2_binette_bins
check_grep "bin_contig_labels: header" "test_sample2_contig_bins.tsv" "^contig	bin$"
check_grep "bin_contig_labels: contig in its bin" "test_sample2_contig_bins.tsv" "^k141_0	binette_bin1$"
check_grep "bin_contig_labels: second bin" "test_sample2_contig_bins.tsv" "^k141_231	binette_bin2$"
check_grep "bin_contig_labels: contig in two bins warned" "../bin_contig_labels.log" "contig k141_0 of bin binette_bin2 is already in bin binette_bin1"
if [ "$(grep -c '^k141_0	' "${BCL_WORK}/test_sample2_contig_bins.tsv")" = "1" ]; then
    echo "✓ bin_contig_labels: a shared contig counts once"; PASS=$((PASS + 1))
else
    echo "✗ bin_contig_labels: a shared contig counts once"; FAIL=$((FAIL + 1))
fi
NEXT_WORKDIR="${BCL_WORK}"
expect_pass "bin_contig_labels: marker-only bins dir" bin_contig_labels.py sampleempty empty_bins
if [ "$(cat "${BCL_WORK}/sampleempty_contig_bins.tsv")" = "contig	bin" ]; then
    echo "✓ bin_contig_labels: marker-only bins dir gives a header-only table"; PASS=$((PASS + 1))
else
    echo "✗ bin_contig_labels: marker-only bins dir gives a header-only table"; FAIL=$((FAIL + 1))
fi
NEXT_WORKDIR="${BCL_WORK}"
expect_fail "bin_contig_labels: missing bins dir argument fails" bin_contig_labels.py test_sample2
ACB_WORK="${WORK}/arg_contig_level_report_bins"
mkdir -p "${ACB_WORK}"
cp "${ACR_WORK}"/test_sample2_* "${BCL_WORK}/test_sample2_contig_bins.tsv" "${BCL_WORK}"/test_sample2_gtdbtk_*.tsv "${ACB_WORK}/"
rm -f "${ACB_WORK}/Contig_tax_and_arg_prediction.tsv"
NEXT_WORKDIR="${ACB_WORK}"
expect_pass "arg_contig_level_report: with bin tables" arg_contig_level_report.py
check_grep "arg_contig_level_report: bin columns after the NCBI ranks" "Contig_tax_and_arg_prediction.tsv" \
    "^Sample	contig_name	contig_length	GC	coverage	superkingdom	phylum	class	order	family	genus	species	taxonomy_rank	taxonomy_score	bin	bin_gtdb_classification	ARG_gene	"
check_grep "arg_contig_level_report: binned contig keeps its Metabuli ranks and gets the GTDB string" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_0	.*	Acinetobacter baumannii	species	0.95811	binette_bin1	d__Bacteria;.*;s__Acinetobacter baumannii	ACRA	"
check_grep "arg_contig_level_report: unbinned contig has empty bin columns" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_456	.*	unclassified	0.0			"
check_grep "arg_contig_level_report: bin without a GTDB-Tk row keeps an empty classification" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_231	.*	genus	0.334152	binette_bin2		"

#
# bin_arg_split.py (Phase 9b, M8): DEEPARG_CONTIGS' hits and the assembly's
# Pyrodigal proteins split by bin (the contig table of bin_contig_labels.py
# above) into ARG_FASTA_FORMATTER's per-bin layout; a bin without hits gets a
# header-only .mapping.ARG; genes of unbinned contigs are dropped.
#
echo "--- bin_arg_split.py ---"
BAS_WORK="${WORK}/bin_arg_split"
mkdir -p "${BAS_WORK}"
cp "${BCL_WORK}/test_sample2_contig_bins.tsv" "${BAS_WORK}/"
printf ">k141_0_1 # 1 # 397 # 1 # ID=k141_0_1\nMKAA\n>k141_0_2 # 1000 # 1150 # 1\nMKCC\n>k141_820_1 # 10 # 300 # 1\nMKDD\n>k141_231_1 # 10 # 300 # -1\nMKEE\n>k141_456_1 # 5 # 200 # 1\nMKFF\n" \
    | gzip > "${BAS_WORK}/test_sample2.faa.gz"
{
    cat "${REPO_DIR}/tests/data/reports/test_sample2_contigs_deep_arg.out.mapping.ARG"
    printf "TETM\t1\t50\tk141_456_1\ttetracycline\tX|FEATURES|tetM|tetracycline|tetM\t0.99\t99.0\t50\t100.0\t1e-20\t1\n"
} > "${BAS_WORK}/test_sample2_contigs_deep_arg.out.mapping.ARG"
NEXT_WORKDIR="${BAS_WORK}"
expect_pass "bin_arg_split: hits and proteins split by bin" bin_arg_split.py test_sample2 test_sample2_contig_bins.tsv \
    test_sample2.faa.gz test_sample2_contigs_deep_arg.out.mapping.ARG
check_file "bin_arg_split: bin1 proteins" "test_sample2_bins_proteins/test_sample2_binette_bin1_proteins.faa"
check_file "bin_arg_split: bin2 proteins" "test_sample2_bins_proteins/test_sample2_binette_bin2_proteins.faa"
check_grep "bin_arg_split: bin1 holds its contigs' proteins with Pyrodigal's headers" \
    "test_sample2_bins_proteins/test_sample2_binette_bin1_proteins.faa" "^>k141_820_1 # 10 # 300 # 1$"
check_grep "bin_arg_split: bin1's ACRA hit with DeepARG's header" \
    "test_sample2_deeparg_results/test_sample2_binette_bin1_deep_arg.out.mapping.ARG" "^ACRA	1	397	k141_0_1	multidrug	"
if [ "$(cat "${BAS_WORK}/test_sample2_deeparg_results/test_sample2_binette_bin2_deep_arg.out.mapping.ARG")" = \
     "$(head -1 "${REPO_DIR}/tests/data/reports/test_sample2_contigs_deep_arg.out.mapping.ARG")" ]; then
    echo "✓ bin_arg_split: a bin without hits gets the header only"; PASS=$((PASS + 1))
else
    echo "✗ bin_arg_split: a bin without hits gets the header only"; FAIL=$((FAIL + 1))
fi
if ! grep -rq "k141_456" "${BAS_WORK}/test_sample2_bins_proteins" "${BAS_WORK}/test_sample2_deeparg_results"; then
    echo "✓ bin_arg_split: genes of unbinned contigs are dropped"; PASS=$((PASS + 1))
else
    echo "✗ bin_arg_split: genes of unbinned contigs are dropped"; FAIL=$((FAIL + 1))
fi
NEXT_WORKDIR="${BAS_WORK}"
expect_fail "bin_arg_split: a file that is not a DeepARG table fails" bin_arg_split.py test_sampleX \
    test_sample2_contig_bins.tsv test_sample2.faa.gz test_sample2_contig_bins.tsv

#
# arg_contig_level_report.py with geNomad (Phase 9b, L6 d): synthetic geNomad
# 1.13.0 summaries (real column layout) for three contigs: k141_820 a plasmid
# (MOBQ), k141_231 a virus, k141_0 a chromosome contig with a provirus at
# 200-900 that overlaps its gene k141_0_1 (1-397) but not k141_0_2
# (1000-1150); gene coordinates from a Pyrodigal-style GFF. Empty geNomad
# files (stub outputs) still give the columns.
#
echo "--- arg_contig_level_report.py (geNomad) ---"
ACG_WORK="${WORK}/arg_contig_level_report_genomad"
mkdir -p "${ACG_WORK}"
cp "${REPO_DIR}/tests/data/reports/test_sample2_Blob_table.test_sample2.blobDB.table.txt" \
   "${REPO_DIR}/tests/data/reports/test_sample2_classifications.tsv" "${ACG_WORK}/"
{
    cat "${REPO_DIR}/tests/data/reports/test_sample2_contigs_deep_arg.out.mapping.ARG"
    printf "TETM\t1\t50\tk141_0_2\ttetracycline\tX|FEATURES|tetM|tetracycline|tetM\t0.99\t99.0\t50\t100.0\t1e-20\t1\n"
    printf "MEFA\t1\t50\tk141_820_1\tmacrolide-lincosamide-streptogramin\tX|FEATURES|mefA|MLS|mefA\t0.98\t98.0\t50\t100.0\t1e-20\t1\n"
    printf "BCRA\t1\t50\tk141_231_1\tbacitracin\tX|FEATURES|bcrA|bacitracin|bcrA\t0.97\t97.0\t50\t100.0\t1e-20\t1\n"
} > "${ACG_WORK}/test_sample2_contigs_deep_arg.out.mapping.ARG"
printf "seq_name\tchromosome_score\tplasmid_score\tvirus_score\nk141_0\t0.9500\t0.0300\t0.0200\nk141_820\t0.0100\t0.9800\t0.0100\nk141_231\t0.0100\t0.0100\t0.9800\n" \
    > "${ACG_WORK}/test_sample2_filtered_contigs_aggregated_classification.tsv"
printf "seq_name\tlength\ttopology\tn_genes\tgenetic_code\tplasmid_score\tfdr\tn_hallmarks\tmarker_enrichment\tconjugation_genes\tamr_genes\nk141_820\t1015\tNo terminal repeats\t2\t11\t0.9800\tNA\t1\t5.0\tMOBQ\tNA\n" \
    > "${ACG_WORK}/test_sample2_filtered_contigs_plasmid_summary.tsv"
printf "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\nk141_231\t1018\tNo terminal repeats\tNA\t2\t11\t0.9800\tNA\t2\t6.4\tViruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;\nk141_0|provirus_200_900\t701\tProvirus\t200-900\t1\t11\t0.9000\tNA\t1\t3.0\tViruses;Duplodnaviria;;;;;\n" \
    > "${ACG_WORK}/test_sample2_filtered_contigs_virus_summary.tsv"
printf "##gff-version  3\nk141_0\tpyrodigal_v3.7.1\tCDS\t1\t397\t50.0\t+\t0\tID=k141_0_1;partial=10;\nk141_0\tpyrodigal_v3.7.1\tCDS\t1000\t1150\t50.0\t+\t0\tID=k141_0_2;partial=00;\nk141_820\tpyrodigal_v3.7.1\tCDS\t10\t300\t50.0\t+\t0\tID=k141_820_1;\nk141_231\tpyrodigal_v3.7.1\tCDS\t10\t300\t50.0\t-\t0\tID=k141_231_1;\n" \
    | gzip > "${ACG_WORK}/test_sample2.gff.gz"
NEXT_WORKDIR="${ACG_WORK}"
expect_pass "arg_contig_level_report: with geNomad tables" arg_contig_level_report.py
check_grep "arg_contig_level_report: geNomad columns before the ARG columns" "Contig_tax_and_arg_prediction.tsv" \
    "	taxonomy_score	mobile_class	plasmid_score	virus_score	provirus	virus_taxonomy	plasmid_conjugation_genes	ARG_gene	"
check_grep "arg_contig_level_report: ARG gene overlapping a provirus" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_0	.*	chromosome	0.0300	0.0200	200-900	Viruses;Duplodnaviria;;;;;		ACRA	"
check_grep "arg_contig_level_report: ARG gene outside the provirus" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_0	.*	chromosome	0.0300	0.0200				TETM	"
check_grep "arg_contig_level_report: plasmid contig with its conjugation genes" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_820	.*	plasmid	0.9800	0.0100			MOBQ	MEFA	"
check_grep "arg_contig_level_report: virus contig with its taxonomy" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_231	.*	virus	0.0100	0.9800		Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;		BCRA	"
check_grep "arg_contig_level_report: contig without a geNomad call is chromosome" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_456	.*	chromosome						"
ACGE_WORK="${WORK}/arg_contig_level_report_genomad_empty"
mkdir -p "${ACGE_WORK}"
cp "${ACR_WORK}"/test_sample2_* "${ACGE_WORK}/"
rm -f "${ACGE_WORK}/Contig_tax_and_arg_prediction.tsv"
touch "${ACGE_WORK}/test_sample2_filtered_contigs_aggregated_classification.tsv" \
      "${ACGE_WORK}/test_sample2_filtered_contigs_plasmid_summary.tsv" \
      "${ACGE_WORK}/test_sample2_filtered_contigs_virus_summary.tsv"
echo "" | gzip > "${ACGE_WORK}/test_sample2.gff.gz"
NEXT_WORKDIR="${ACGE_WORK}"
expect_pass "arg_contig_level_report: empty geNomad files (stub) still give the columns" arg_contig_level_report.py
check_grep "arg_contig_level_report: empty geNomad files give chromosome" "Contig_tax_and_arg_prediction.tsv" \
    "^test_sample2	k141_0	.*	chromosome						ACRA	"

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
