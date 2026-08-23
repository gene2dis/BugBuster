#!/bin/bash
#
# Regression tests for audit item #4:
#   - bin/taxonomy_report.py must parse the raw headerless 6-column kraken2
#     reports the nf-core KRAKEN2 module actually produces, and
#   - bin/taxonomy_phyloseq.py must parse the kraken-style bracken reports
#     (BRACKEN.out.txt) fed to it, merging samples on the union of taxa with
#     zeros (not NaN) and a taxonomy table covering every sample.
#
# The scripts run inside the same pinned image the report processes use,
# against the committed fixtures in tests/data/taxonomy/.
#
# Requirements: docker (present on GitHub ubuntu-latest runners).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
FIXTURES="${REPO_DIR}/tests/data/taxonomy"

# Use the exact image pinned in the module so the test cannot drift from it
REPORT_IMG=$(grep -o "'[^']*community.wave.seqera.io/library/python_pandas[^']*'" "${REPO_DIR}/modules/local/taxonomy_report/main.nf" \
    | tr -d "'" | grep -v '^oras://' | head -1)

if [ -z "${REPORT_IMG}" ]; then
    echo "ERROR: could not extract pinned container image from modules/local/taxonomy_report/main.nf"
    exit 1
fi

echo "=== Taxonomy script test suite ==="
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

run_script() {
    # run_script <workdir> <script-name> [args...] — inside the report image
    local workdir="$1" script="$2"
    shift 2
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e MPLCONFIGDIR=/tmp \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" \
        -v "${workdir}:/taxwork" -w /taxwork \
        "${REPORT_IMG}" python3 "/pipeline_bin/${script}" "$@"
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
# taxonomy_report.py — kraken2 branch (audit #4)
#
echo "--- taxonomy_report.py (kraken2) ---"
KR_WORK="${WORK}/taxonomy_report_kraken2"
mkdir -p "${KR_WORK}"
cp "${FIXTURES}/sample1.kraken2.report.txt" "${FIXTURES}/sample2.kraken2.report.txt" "${KR_WORK}/"
cp "${FIXTURES}/Reads_report_fixture.csv" "${KR_WORK}/input_reads_report.csv"
LAST_WORKDIR="${KR_WORK}"
if run_script "${KR_WORK}" taxonomy_report.py --profiler kraken2 \
    --reports sample1.kraken2.report.txt sample2.kraken2.report.txt \
    --reads-report input_reads_report.csv --db-name standard-8 --output-dir . \
    > "${KR_WORK}.log" 2>&1; then
    echo "✓ kraken2 report: raw nf-core reports parsed"
    PASS=$((PASS + 1))
else
    echo "✗ kraken2 report: raw nf-core reports parsed — failed:"
    tail -5 "${KR_WORK}.log" | sed 's/^/    /'
    FAIL=$((FAIL + 1))
fi
check_grep "kraken2 report: merged header keeps QC columns" "Reads_report.csv" \
    "Id,Raw reads,Fastp,Clean reads (contaminants removed),Kraken DB,Unclassified,Classified"
check_grep "kraken2 report: sample1 classified fractions" "Reads_report.csv" "^sample1,.*,standard-8,12.5,87.5$"
check_grep "kraken2 report: sample2 (no U row) fully classified" "Reads_report.csv" "^sample2,.*,standard-8,0.0,100.0$"
check_file "kraken2 report: classification plot created" "Kraken_plot.png"

# Regression (symlink-corruption fix): the staged input reads report must
# not be modified by the run
if cmp -s "${KR_WORK}/input_reads_report.csv" "${FIXTURES}/Reads_report_fixture.csv"; then
    echo "✓ kraken2 report: input reads report untouched"
    PASS=$((PASS + 1))
else
    echo "✗ kraken2 report: input reads report was modified by the script"
    FAIL=$((FAIL + 1))
fi

# Garbled report (not 6 tab-separated columns) must fail loudly
GARBLED="${WORK}/garbled"
mkdir -p "${GARBLED}"
printf 'just\tthree\tcolumns\n' > "${GARBLED}/bad.kraken2.report.txt"
cp "${FIXTURES}/Reads_report_fixture.csv" "${GARBLED}/input_reads_report.csv"
NEXT_WORKDIR="${GARBLED}"
expect_fail "kraken2 report: garbled report fails" taxonomy_report.py --profiler kraken2 \
    --reports bad.kraken2.report.txt \
    --reads-report input_reads_report.csv --db-name standard-8 --output-dir .

#
# taxonomy_report.py — sourmash branch non-regression (format unchanged)
#
echo "--- taxonomy_report.py (sourmash) ---"
SM_WORK="${WORK}/taxonomy_report_sourmash"
mkdir -p "${SM_WORK}"
cp "${FIXTURES}/sample1_species_gtdb_report.tsv" "${SM_WORK}/"
cp "${FIXTURES}/Reads_report_fixture.csv" "${SM_WORK}/input_reads_report.csv"
NEXT_WORKDIR="${SM_WORK}"
expect_pass "sourmash report: headered TSV still parsed" taxonomy_report.py --profiler sourmash \
    --reports sample1_species_gtdb_report.tsv \
    --reads-report input_reads_report.csv --db-name gtdb --output-dir .
check_grep "sourmash report: classified fraction merged" "Reads_report.csv" "^sample1,.*,gtdb,0.25,0.75$"
check_file "sourmash report: classification plot created" "sourmash_tax_classified_reads.png"

#
# taxonomy_phyloseq.py — kraken2 branch (audit #4; merge bugs were #12)
#
echo "--- taxonomy_phyloseq.py (kraken2) ---"
PH_WORK="${WORK}/taxonomy_phyloseq_kraken2"
mkdir -p "${PH_WORK}"
cp "${FIXTURES}/sample1.kraken2.report_bracken.txt" "${FIXTURES}/sample2.kraken2.report_bracken.txt" "${PH_WORK}/"
NEXT_WORKDIR="${PH_WORK}"
expect_pass "kraken2 phyloseq: bracken kraken-style reports parsed" taxonomy_phyloseq.py --profiler kraken2 \
    --input-files sample1.kraken2.report_bracken.txt sample2.kraken2.report_bracken.txt \
    --db-name standard-8 --output-dir . --format both \
    --plot-levels Phylum,Family,Genus,Species --top-n 10
check_grep "kraken2 phyloseq: both sample columns in OTU table" "kraken2_standard-8_otu_table.tsv" \
    "sample1	sample2"
check_grep "kraken2 phyloseq: shared species counted in both samples" "kraken2_standard-8_otu_table.tsv" \
    "^562	610.0	380.0$"
check_grep "kraken2 phyloseq: sample1-only species is 0.0 (not NaN) in sample2" "kraken2_standard-8_otu_table.tsv" \
    "^1423	265.0	0.0$"
check_grep "kraken2 phyloseq: sample2-only species is 0.0 (not NaN) in sample1" "kraken2_standard-8_otu_table.tsv" \
    "^1280	0.0	420.0$"
check_grep "kraken2 phyloseq: tax table has full lineage for sample1 taxon" "kraken2_standard-8_tax_table.tsv" \
    "^1423	Bacteria	Bacillota	Bacilli	Bacillales	Bacillaceae	Bacillus	Bacillus subtilis$"
check_grep "kraken2 phyloseq: tax table covers sample2-only taxon (union)" "kraken2_standard-8_tax_table.tsv" \
    "^1280	Bacteria	Bacillota	Bacilli	Bacillales	Staphylococcaceae	Staphylococcus	Staphylococcus aureus$"
check_file "kraken2 phyloseq: sample metadata written" "kraken2_standard-8_sample_metadata.tsv"
check_file "kraken2 phyloseq: HDF5 written" "kraken2_standard-8_phyloseq_data.h5"
check_file "kraken2 phyloseq: species plot in output root" "kraken2_standard-8_species_bar.png"
check_file "kraken2 phyloseq: phylum plot in output root" "kraken2_standard-8_phylum_bar.png"

# Garbled input must fail loudly
GARBLED_PH="${WORK}/garbled_phyloseq"
mkdir -p "${GARBLED_PH}"
printf 'not\ta\tkraken\treport\n' > "${GARBLED_PH}/bad.kraken2.report_bracken.txt"
NEXT_WORKDIR="${GARBLED_PH}"
expect_fail "kraken2 phyloseq: garbled report fails" taxonomy_phyloseq.py --profiler kraken2 \
    --input-files bad.kraken2.report_bracken.txt \
    --db-name standard-8 --output-dir . --format tables \
    --plot-levels Species --top-n 10

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
