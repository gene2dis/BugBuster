#!/bin/bash
#
# Regression tests for audit items #4 and #12:
#   - bin/taxonomy_report.py must parse the raw headerless 6-column kraken2
#     reports the nf-core KRAKEN2 module actually produces (#4),
#   - bin/taxonomy_phyloseq.py must parse the kraken-style bracken reports
#     (BRACKEN.out.txt) fed to it, merging samples on the union of taxa with
#     zeros (not NaN) and a taxonomy table covering every sample (#4), and
#   - the sourmash branch must strip module filename decorations from sample
#     ids, skip a no-match sample (the SOURMASH module's header-only fallback
#     CSV), and fail loudly when no sample has any match (#12).
#   - the sylph branch must merge only the four standard report columns into
#     Reads_report.csv, build the relative-abundance OTU table from the t__
#     genome rows of the .sylphmpa profiles, the estimated read pairs table
#     from sequence_abundance x the report's Read pairs, skip a no-hit sample
#     (SYLPH_TAXPROF's header-only fallback) and fail loudly otherwise.
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
# REPORT_TEST_IMG=<image>: run the scripts in that image instead of the pinned one (to try an image before pinning it)
REPORT_IMG=${REPORT_TEST_IMG:-${REPORT_IMG}}

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

#
# taxonomy_phyloseq.py — sourmash branch (audit #12)
#
echo "--- taxonomy_phyloseq.py (sourmash) ---"

# The no-match fallback CSV is generated from the module itself so this test
# cannot drift from what the pipeline really produces
FALLBACK_HEADER=$(grep -o '"intersect_bp,[^"]*"' "${REPO_DIR}/modules/local/sourmash/main.nf" | tr -d '"' | head -1)
if [ -z "${FALLBACK_HEADER}" ]; then
    echo "✗ sourmash phyloseq: could not extract fallback CSV header from modules/local/sourmash/main.nf"
    FAIL=$((FAIL + 1))
else
    echo "✓ sourmash phyloseq: fallback CSV header extracted from module"
    PASS=$((PASS + 1))
fi

# The committed header-only fixture (used by the nf-test suite) must match
# the module's fallback header
if [ "$(head -1 "${FIXTURES}/sampleB_smgather_gtdb.with-lineages.csv")" = "${FALLBACK_HEADER}" ]; then
    echo "✓ sourmash phyloseq: committed no-match fixture matches module fallback"
    PASS=$((PASS + 1))
else
    echo "✗ sourmash phyloseq: sampleB fixture header drifted from module fallback"
    FAIL=$((FAIL + 1))
fi

SMPH_WORK="${WORK}/taxonomy_phyloseq_sourmash"
mkdir -p "${SMPH_WORK}"
cp "${FIXTURES}/sampleA_smgather_gtdb.with-lineages.csv" "${SMPH_WORK}/"
echo "${FALLBACK_HEADER}" > "${SMPH_WORK}/sampleB_smgather_gtdb.with-lineages.csv"
NEXT_WORKDIR="${SMPH_WORK}"
expect_pass "sourmash phyloseq: match + no-match samples parsed" taxonomy_phyloseq.py --profiler sourmash \
    --input-files sampleA_smgather_gtdb.with-lineages.csv sampleB_smgather_gtdb.with-lineages.csv \
    --db-name gtdb --output-dir . --format tables \
    --plot-levels Phylum,Species --top-n 10
check_grep "sourmash phyloseq: sample id has no .with-lineages suffix" "sourmash_gtdb_otu_table.tsv" \
    "^name	sampleA$"
check_grep "sourmash phyloseq: abundance computed from scaled column" "sourmash_gtdb_otu_table.tsv" \
    "^GCA_000000001.1	200\(\.0\)\{0,1\}$"
check_grep "sourmash phyloseq: lineage split into ranks" "sourmash_gtdb_tax_table.tsv" \
    "^GCA_000000002.1	Bacteria	Bacillota	Bacilli	Bacillales	Bacillaceae	Bacillus	Bacillus subtilis$"
check_grep "sourmash phyloseq: no-match sample warned and skipped" "../taxonomy_phyloseq_sourmash.log" \
    "no sourmash matches for sample sampleB"

# All samples empty (the audit-requested no-matches unit test): loud failure
SMPH_EMPTY="${WORK}/taxonomy_phyloseq_sourmash_empty"
mkdir -p "${SMPH_EMPTY}"
echo "${FALLBACK_HEADER}" > "${SMPH_EMPTY}/sampleB_smgather_gtdb.with-lineages.csv"
NEXT_WORKDIR="${SMPH_EMPTY}"
expect_fail "sourmash phyloseq: all-samples-no-match fails" taxonomy_phyloseq.py --profiler sourmash \
    --input-files sampleB_smgather_gtdb.with-lineages.csv \
    --db-name gtdb --output-dir . --format tables \
    --plot-levels Species --top-n 10
check_grep "sourmash phyloseq: no-match failure message is actionable" "../taxonomy_phyloseq_sourmash_empty.log" \
    "No sourmash matches in any sample"

# A gather CSV missing expected columns (e.g. the old fallback header with
# match_name and no scaled) must fail loudly, not be silently skipped
SMPH_BADHDR="${WORK}/taxonomy_phyloseq_sourmash_badhdr"
mkdir -p "${SMPH_BADHDR}"
printf 'query_filename,match_name,average_abund,lineage\nx.fq,GCA_1 Test,1.0,d__Bacteria\n' \
    > "${SMPH_BADHDR}/sampleC_smgather_gtdb.with-lineages.csv"
NEXT_WORKDIR="${SMPH_BADHDR}"
expect_fail "sourmash phyloseq: malformed gather CSV fails" taxonomy_phyloseq.py --profiler sourmash \
    --input-files sampleC_smgather_gtdb.with-lineages.csv \
    --db-name gtdb --output-dir . --format tables \
    --plot-levels Species --top-n 10

#
# taxonomy_report.py — sylph branch
#
echo "--- taxonomy_report.py (sylph) ---"
SY_WORK="${WORK}/taxonomy_report_sylph"
mkdir -p "${SY_WORK}"
cp "${FIXTURES}/sample1_sylph_gtdb_report.tsv" "${SY_WORK}/"
cp "${FIXTURES}/Reads_report_fixture.csv" "${SY_WORK}/input_reads_report.csv"
NEXT_WORKDIR="${SY_WORK}"
expect_pass "sylph report: SYLPH_PROFILE report parsed" taxonomy_report.py --profiler sylph \
    --reports sample1_sylph_gtdb_report.tsv \
    --reads-report input_reads_report.csv --db-name gtdb --output-dir .
check_grep "sylph report: Sylph DB column in header" "Reads_report.csv" ",Sylph DB,Unclassified,Classified$"
check_grep "sylph report: classified fraction merged, Read pairs dropped" "Reads_report.csv" "^sample1,.*,gtdb,0.2404,0.7596$"
check_file "sylph report: classification plot created" "sylph_tax_classified_reads.png"

# A per-sample report missing a standard column must fail, not leave the
# sample's classification blank in Reads_report.csv
SY_BADREP="${WORK}/taxonomy_report_sylph_badcols"
mkdir -p "${SY_BADREP}"
printf 'Id\tSylph DB\tClassified\nsample1\tgtdb\t0.75\n' > "${SY_BADREP}/sample1_sylph_gtdb_report.tsv"
cp "${FIXTURES}/Reads_report_fixture.csv" "${SY_BADREP}/input_reads_report.csv"
NEXT_WORKDIR="${SY_BADREP}"
expect_fail "sylph report: report missing a column fails" taxonomy_report.py --profiler sylph \
    --reports sample1_sylph_gtdb_report.tsv \
    --reads-report input_reads_report.csv --db-name gtdb --output-dir .
check_grep "sylph report: missing-column failure names the column" "../taxonomy_report_sylph_badcols.log" \
    "missing report columns \['Unclassified'\]"

#
# taxonomy_phyloseq.py — sylph branch
#
echo "--- taxonomy_phyloseq.py (sylph) ---"

# The no-hit fallback header comes from the module, so the committed
# header-only fixture cannot drift from what SYLPH_TAXPROF writes
SY_HEADER=$(grep -o "def header = '[^']*'" "${REPO_DIR}/modules/local/sylph_taxprof/main.nf" \
    | sed "s/^def header = '//; s/'\$//; s/\\\\\\\\t/\t/g")
if [ "$(sed -n 2p "${FIXTURES}/sampleB_sylph_gtdb.sylphmpa")" = "${SY_HEADER}" ] \
    && [ "$(sed -n 2p "${FIXTURES}/sampleA_sylph_gtdb.sylphmpa")" = "${SY_HEADER}" ]; then
    echo "✓ sylph phyloseq: fixture headers match the SYLPH_TAXPROF fallback header"
    PASS=$((PASS + 1))
else
    echo "✗ sylph phyloseq: .sylphmpa fixture header drifted from modules/local/sylph_taxprof/main.nf"
    FAIL=$((FAIL + 1))
fi

SYPH_WORK="${WORK}/taxonomy_phyloseq_sylph"
mkdir -p "${SYPH_WORK}"
cp "${FIXTURES}"/sample[AB]_sylph_gtdb.sylphmpa "${FIXTURES}"/sample[AB]_sylph_gtdb_report.tsv "${SYPH_WORK}/"
NEXT_WORKDIR="${SYPH_WORK}"
expect_pass "sylph phyloseq: hit + no-hit samples parsed" taxonomy_phyloseq.py --profiler sylph \
    --input-files sampleA_sylph_gtdb.sylphmpa sampleA_sylph_gtdb_report.tsv \
        sampleB_sylph_gtdb.sylphmpa sampleB_sylph_gtdb_report.tsv \
    --db-name gtdb --output-dir . --format tables \
    --plot-levels Phylum,Species --top-n 10
check_grep "sylph phyloseq: OTU table keyed by t__ accession, sample id from file name" "sylph_gtdb_otu_table.tsv" \
    "^taxon_id	sampleA$"
check_grep "sylph phyloseq: OTU table holds relative_abundance" "sylph_gtdb_otu_table.tsv" \
    "^GCF_000009045.1	22.3376$"
check_grep "sylph phyloseq: only genome rows (3 genomes + header)" "sylph_gtdb_otu_table.tsv" \
    "^GCA_047199095.1	37.8619$"
if [ "$(wc -l < "${SYPH_WORK}/sylph_gtdb_otu_table.tsv" 2>/dev/null)" = "4" ]; then
    echo "✓ sylph phyloseq: higher-rank rows not counted as taxa"
    PASS=$((PASS + 1))
else
    echo "✗ sylph phyloseq: OTU table should have 3 genome rows"
    FAIL=$((FAIL + 1))
fi
check_grep "sylph phyloseq: read pairs = sequence_abundance x Read pairs" "sylph_gtdb_read_counts.tsv" \
    "^GCA_047199115.1	24465$"
check_grep "sylph phyloseq: read pairs for B. subtilis" "sylph_gtdb_read_counts.tsv" \
    "^GCF_000009045.1	12507$"
check_grep "sylph phyloseq: lineage split by rank prefix" "sylph_gtdb_tax_table.tsv" \
    "^GCF_000009045.1	Bacteria	Bacillota	Bacilli	Bacillales	Bacillaceae	Bacillus	Bacillus subtilis$"
check_grep "sylph phyloseq: no-hit sample warned and skipped" "../taxonomy_phyloseq_sylph.log" \
    "no sylph genomes for sample sampleB"

# A custom database label (its file name) may itself contain 'sylph': the
# sample id is the file name minus the exact _sylph_<db>.sylphmpa suffix
SYPH_LABEL="${WORK}/taxonomy_phyloseq_sylph_label"
mkdir -p "${SYPH_LABEL}"
cp "${FIXTURES}/sampleA_sylph_gtdb.sylphmpa" "${SYPH_LABEL}/sampleA_sylph_my_sylph_db.sylphmpa"
cp "${FIXTURES}/sampleA_sylph_gtdb_report.tsv" "${SYPH_LABEL}/sampleA_sylph_my_sylph_db_report.tsv"
NEXT_WORKDIR="${SYPH_LABEL}"
expect_pass "sylph phyloseq: database label containing 'sylph' parsed" taxonomy_phyloseq.py --profiler sylph \
    --input-files sampleA_sylph_my_sylph_db.sylphmpa sampleA_sylph_my_sylph_db_report.tsv \
    --db-name my_sylph_db --output-dir . --format tables \
    --plot-levels Species --top-n 10
check_grep "sylph phyloseq: sample id is the name minus the exact suffix" "sylph_my_sylph_db_otu_table.tsv" \
    "^taxon_id	sampleA$"

# A profile not named for the --db-name database must fail, not be mislabelled
NEXT_WORKDIR="${SYPH_LABEL}"
expect_fail "sylph phyloseq: profile of another database label fails" taxonomy_phyloseq.py --profiler sylph \
    --input-files sampleA_sylph_my_sylph_db.sylphmpa sampleA_sylph_my_sylph_db_report.tsv \
    --db-name gtdb --output-dir . --format tables \
    --plot-levels Species --top-n 10

# A non-integer Read pairs value (e.g. awk's 1.23457e+07) must fail, not be
# truncated into a wrong count
SYPH_BADPAIRS="${WORK}/taxonomy_phyloseq_sylph_badpairs"
mkdir -p "${SYPH_BADPAIRS}"
cp "${FIXTURES}/sampleA_sylph_gtdb.sylphmpa" "${SYPH_BADPAIRS}/"
printf 'Id\tSylph DB\tUnclassified\tClassified\tRead pairs\nsampleA\tgtdb\t0.000000\t1.000000\t1.23457e+07\n' \
    > "${SYPH_BADPAIRS}/sampleA_sylph_gtdb_report.tsv"
NEXT_WORKDIR="${SYPH_BADPAIRS}"
expect_fail "sylph phyloseq: non-integer Read pairs fails" taxonomy_phyloseq.py --profiler sylph \
    --input-files sampleA_sylph_gtdb.sylphmpa sampleA_sylph_gtdb_report.tsv \
    --db-name gtdb --output-dir . --format tables \
    --plot-levels Species --top-n 10

# All samples without hits: loud failure
SYPH_EMPTY="${WORK}/taxonomy_phyloseq_sylph_empty"
mkdir -p "${SYPH_EMPTY}"
cp "${FIXTURES}/sampleB_sylph_gtdb.sylphmpa" "${FIXTURES}/sampleB_sylph_gtdb_report.tsv" "${SYPH_EMPTY}/"
NEXT_WORKDIR="${SYPH_EMPTY}"
expect_fail "sylph phyloseq: all-samples-no-hit fails" taxonomy_phyloseq.py --profiler sylph \
    --input-files sampleB_sylph_gtdb.sylphmpa sampleB_sylph_gtdb_report.tsv \
    --db-name gtdb --output-dir . --format tables \
    --plot-levels Species --top-n 10
check_grep "sylph phyloseq: no-hit failure message is actionable" "../taxonomy_phyloseq_sylph_empty.log" \
    "No sylph genomes in any sample"

# A profile with hits but no matching report (no Read pairs) must fail
SYPH_NOREP="${WORK}/taxonomy_phyloseq_sylph_noreport"
mkdir -p "${SYPH_NOREP}"
cp "${FIXTURES}/sampleA_sylph_gtdb.sylphmpa" "${SYPH_NOREP}/"
NEXT_WORKDIR="${SYPH_NOREP}"
expect_fail "sylph phyloseq: profile without its report fails" taxonomy_phyloseq.py --profiler sylph \
    --input-files sampleA_sylph_gtdb.sylphmpa \
    --db-name gtdb --output-dir . --format tables \
    --plot-levels Species --top-n 10

# Not a sylph-tax profile (no #SampleID line) must fail
SYPH_BAD="${WORK}/taxonomy_phyloseq_sylph_bad"
mkdir -p "${SYPH_BAD}"
printf 'clade_name\trelative_abundance\nd__Bacteria\t100\n' > "${SYPH_BAD}/sampleC_sylph_gtdb.sylphmpa"
cp "${FIXTURES}/sampleA_sylph_gtdb_report.tsv" "${SYPH_BAD}/sampleC_sylph_gtdb_report.tsv"
NEXT_WORKDIR="${SYPH_BAD}"
expect_fail "sylph phyloseq: malformed profile fails" taxonomy_phyloseq.py --profiler sylph \
    --input-files sampleC_sylph_gtdb.sylphmpa sampleC_sylph_gtdb_report.tsv \
    --db-name gtdb --output-dir . --format tables \
    --plot-levels Species --top-n 10

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
