#!/bin/bash
#
# Tests for bin/binette_quality_reports.py: Binette's CheckM2 reports
# (input_bins_quality_reports/*.tsv, final_bins_quality_reports.tsv) written in
# CHECKM2_BATCH's per-binner format, which BIN_QUALITY_REPORT and BIN_SUMMARY
# read (performance review M7). Synthetic fixture built here; no containers.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
SCRIPT="${REPO_DIR}/bin/binette_quality_reports.py"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0
ok()   { echo "✓ $1"; PASS=$((PASS + 1)); }
bad()  { echo "✗ $1"; FAIL=$((FAIL + 1)); }
check() { if eval "$2"; then ok "$1"; else bad "$1"; fi; }

HEADER=$'Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tTotal_Contigs\tMax_Contig_Length\tAdditional_Notes'
BINETTE_HEADER=$'name\torigin\tis_original\toriginal_name\tcompleteness\tcontamination\tscore\tcheckm2_model\tsize\tN50\tcoding_density\tcontig_count'

# Fixture: two binners with bins, COMEBin with only a FAILED.txt note (staged
# under its nested basename), and Binette's output for the two sets
make_fixture() {
    local d="$1"
    mkdir -p "$d/S1_metabat_bins" "$d/S1_semibin_output_bins" "$d/comebin_res_bins" \
        "$d/binette_out/input_bins_quality_reports" "$d/binette_out/temporary_files" "$d/final"
    # MetaBAT2 bins: c1 (GC 50 %, 8 bp) + c2 (GC 100 %, 4 bp) | c3
    printf '>c1\nACGTACGT\n>c2\nGGCC\n' > "$d/S1_metabat_bins/S1_metabat_bin.1.fa"
    printf '>c3\nAAAA\n' > "$d/S1_metabat_bins/S1_metabat_bin.2.fa"
    printf '>c1\nACGTACGT\n' > "$d/S1_semibin_output_bins/SemiBin_1.fa"
    echo "COMEBIN failed" > "$d/comebin_res_bins/FAILED.txt"
    printf '>c1\nACGTACGT\n>c2\nGGCC\n' > "$d/final/binette_bin1.fa"
    {
        echo "${BINETTE_HEADER}"
        printf 'S1_metabat_bin.2\tS1_metabat_bins\tTrue\tS1_metabat_bin.2\t10.0\t1.0\t8.0\tGradient Boost (General Model)\t4\t4\t0.5\t1\n'
        printf 'S1_metabat_bin.1\tS1_metabat_bins\tTrue\tS1_metabat_bin.1\t95.5\t2.1\t91.3\tNeural Network (Specific Model)\t12\t8\t0.9\t2\n'
    } > "$d/binette_out/input_bins_quality_reports/input_bins_1.S1_metabat_bins.tsv"
    {
        echo "${BINETTE_HEADER}"
        printf 'SemiBin_1\tS1_semibin_output_bins\tTrue\tSemiBin_1\t60.0\t3.0\t54.0\tNeural Network (Specific Model)\t8\t8\t0.8\t1\n'
    } > "$d/binette_out/input_bins_quality_reports/input_bins_2.S1_semibin_output_bins.tsv"
    {
        echo "${BINETTE_HEADER}"
        printf 'binette_bin1\tS1_metabat_bins\tTrue\tS1_metabat_bin.1\t95.5\t2.1\t91.3\tNeural Network (Specific Model)\t12\t8\t0.9\t2\n'
    } > "$d/binette_out/final_bins_quality_reports.tsv"
    # Prodigal-format proteins: two genes on c1 (3 and 5 aa), one on c2 (4 aa), none on c3
    printf '>c1_1 # 1 # 9 # 1\nMKV*\n>c1_2 # 10 # 24 # 1\nMKVLA\n>c2_1 # 1 # 12 # 1\nMKVL\n' \
        | gzip -c > "$d/binette_out/temporary_files/assembly_proteins.faa.gz"
}

run() {
    # run <dir> <args...>: runs in <dir>, stderr to <dir>/stderr.log
    local d="$1"; shift
    (cd "$d" && python3 "${SCRIPT}" "$@" 2> stderr.log)
}

echo "=== binette_quality_reports.py test suite ==="

D="${TMP_DIR}/happy"
make_fixture "$D"
if run "$D" --sample S1 --binette-out binette_out --final-bins final \
        S1_metabat_bins S1_semibin_output_bins comebin_res_bins; then
    ok "happy path exits 0"
else
    bad "happy path exits 0"; sed 's/^/    /' "$D/stderr.log"
fi
for b in metabat semibin comebin binette; do
    check "report for ${b} written with CheckM2 1.1.0's header" \
        "[ \"\$(head -1 '$D/S1_${b}_quality_report.tsv')\" = \"\${HEADER}\" ]"
done
check "binner without bins (COMEBin FAILED.txt) gets a header-only report" \
    "[ \$(wc -l < '$D/S1_comebin_quality_report.tsv') -eq 1 ]"
check "rows sorted by bin name" \
    "[ \"\$(cut -f1 '$D/S1_metabat_quality_report.tsv' | tail -n +2 | tr '\n' ' ')\" = 'S1_metabat_bin.1 S1_metabat_bin.2 ' ]"
ROW=$(grep '^S1_metabat_bin.1' "$D/S1_metabat_quality_report.tsv")
EXPECTED=$'S1_metabat_bin.1\t95.5\t2.1\tNeural Network (Specific Model)\t\t0.9\t8\t4.0\t12\t0.67\t3\t2\t8\t'
check "Binette columns carried over, gene counts/length from its proteins, GC and longest contig from the FASTA" \
    "[ \"\${ROW}\" = \"\${EXPECTED}\" ]"
check "bin with no gene gets 0 genes and length 0" \
    "grep -q \$'^S1_metabat_bin.2\t10.0\t1.0\t.*\t0\t4\t0.00\t0\t1\t4\t\$' '$D/S1_metabat_quality_report.tsv'"
check "refined bins reported as <sample>_binette_quality_report.tsv" \
    "grep -q \$'^binette_bin1\t95.5\t2.1' '$D/S1_binette_quality_report.tsv'"

D="${TMP_DIR}/not_run"
make_fixture "$D"
check "Binette not run (no --binette-out): exits 0" \
    "run '$D' --sample S1 --final-bins final S1_metabat_bins S1_semibin_output_bins"
check "Binette not run: every report header-only" \
    "[ \$(cat '$D/S1_metabat_quality_report.tsv' '$D/S1_semibin_quality_report.tsv' '$D/S1_binette_quality_report.tsv' | wc -l) -eq 3 ]"

D="${TMP_DIR}/no_proteins"
make_fixture "$D"
rm "$D/binette_out/temporary_files/assembly_proteins.faa.gz"
check "missing protein file: exits 0" \
    "run '$D' --sample S1 --binette-out binette_out --final-bins final S1_metabat_bins S1_semibin_output_bins"
check "missing protein file: warns" "grep -q 'gene counts left empty' '$D/stderr.log'"
check "missing protein file: gene columns empty" \
    "grep -q \$'^S1_metabat_bin.1\t95.5\t2.1\t.*\t0.9\t8\t\t12\t0.67\t\t2\t8\t\$' '$D/S1_metabat_quality_report.tsv'"

D="${TMP_DIR}/missing_bin"
make_fixture "$D"
rm "$D/S1_metabat_bins/S1_metabat_bin.2.fa"
if run "$D" --sample S1 --binette-out binette_out --final-bins final S1_metabat_bins S1_semibin_output_bins; then
    bad "bin reported by Binette but absent from its directory fails"
else
    check "bin reported by Binette but absent from its directory fails" "grep -q 'not in S1_metabat_bins' '$D/stderr.log'"
fi

D="${TMP_DIR}/unknown_dir"
make_fixture "$D"
mv "$D/S1_semibin_output_bins" "$D/S1_random_bins"
if run "$D" --sample S1 --binette-out binette_out --final-bins final S1_metabat_bins S1_random_bins; then
    bad "unrecognized binner directory fails"
else
    check "unrecognized binner directory fails" "grep -q 'unrecognized binner directory' '$D/stderr.log'"
fi

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
