#!/bin/bash
#
# Direct tests for bin/aggregate_functions.py (functional annotation
# T4+T5+T6, design doc Sections 4.3, 4.8, 5, 6, 10.2). Asserts on the
# committed fixtures in tests/data/functional/ (regenerate with
# make_functional_fixtures.sh):
#   - assembly mode: TPM sums to 1e6 per sample; every predicted gene appears
#     in gene_abundance including count-0 genes; a gene with two KOs
#     contributes its full abundance to each (intentional double-counting);
#     multi-letter COG splits per character; '-' fields yield no rows;
#     unannotated genes are absent from 5.1 but present in 5.2; annotated
#     fractions by count and abundance are exact
#   - CPGE (T5): hand-checked cpge = RPK/GE values for the sample with an
#     ags table; empty cpge fields, blank wide-matrix columns and an
#     'unavailable' ags_and_ge.tsv row for samples without one (the TPM-only
#     fallback); running without --ags at all still succeeds
#   - coassembly mode: one shared annotations/GFF, per-sample counts;
#     5.1 sample_id is 'coassembly'; zero-total sample gets tpm 0 and an
#     empty fraction_by_abundance; per-sample ags joins as under assembly
#   - empty gene set: header-only inputs aggregate without crashing
#   - the version-aware layout guard fails loudly (10.2 mandates testing
#     this deliberately): renamed / 21-column / 23-column headers, unknown
#     or missing emapper version, a disagreeing '## emapper-' line, and
#     mis-paired gene id sets all exit non-zero
#   - the ags guard fails loudly on a wrong header, an unknown sample id, a
#     non-positive genome_equivalents, and a filename/embedded-id mismatch
#   - Pfam (Q15): '<name>_<start>_<end>' values stripped to the name, a
#     repeated domain counted once per gene; a suffix-less value fails
#   - dbCAN (T6): db=dbcan rows with per-row run_dbcan provenance in 5.1,
#     backend=run_dbcan rows in 5.3 mirroring the gene's TPM/CPGE, the
#     'recommended' vs 'any' consensus policies, separate cazy vs cazy_dbcan
#     wide matrices (never merged), cazy fraction counting either backend
#     plus a cazy_dbcan row; without --dbcan the cazy_dbcan outputs are
#     header-only/zero. The dbCAN guard fails loudly on a corrupt overview
#     header, a bad row, an unknown run_dbcan version, --dbcan without its
#     versions.yml, a stray gene id, a mismatched sample set, and a
#     duplicate overview
#
# The script runs inside the same pinned image the module uses.
# Requirements: docker (present on GitHub ubuntu-latest runners).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
FIXTURES="${REPO_DIR}/tests/data/functional"
GFF_FIXTURE="${REPO_DIR}/tests/data/gff/test_contig1_genes.gff.gz"
EMPTY_GFF_FIXTURE="${REPO_DIR}/tests/data/gff/test_empty_genes.gff.gz"

# Use the exact image pinned in the module so the test cannot drift from it
REPORT_IMG=$(grep -o "'[^']*community.wave.seqera.io/library/python_pandas[^']*'" "${REPO_DIR}/modules/local/aggregate_functions/main.nf" \
    | tr -d "'" | grep -v '^oras://' | head -1)

if [ -z "${REPORT_IMG}" ]; then
    echo "ERROR: could not extract pinned container image from modules/local/aggregate_functions/main.nf"
    exit 1
fi

echo "=== aggregate_functions.py test suite ==="
echo "Image: ${REPORT_IMG}"
echo ""

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT
WORK="${TMP_DIR}/work"
mkdir -p "${WORK}"

PASS=0
FAIL=0
LAST_WORKDIR=""

run_script() {
    # run_script <workdir> [script args...] — inside the pinned image
    local workdir="$1"
    shift
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e MPLCONFIGDIR=/tmp \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" \
        -v "${workdir}:/aggwork" -w /aggwork \
        "${REPORT_IMG}" python3 /pipeline_bin/aggregate_functions.py "$@"
}

expect_pass() {
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
        echo "✗ ${desc} — expected failure, got success"
        FAIL=$((FAIL + 1))
    else
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    fi
}

check() {
    local desc="$1" cond="$2"
    if eval "${cond}"; then
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${desc}"
        FAIL=$((FAIL + 1))
    fi
}

check_grep() {
    # check_grep <desc> <file rel to LAST_WORKDIR> <pattern>
    check "$1" "grep -q -- '$3' '${LAST_WORKDIR}/$2'"
}

seed_assembly_inputs() {
    # Two-sample assembly-mode input set in $1
    local dir="$1"
    mkdir -p "${dir}"
    cp "${FIXTURES}/sampleA.featureCounts.txt" "${FIXTURES}/sampleB.featureCounts.txt" "${dir}/"
    cp "${FIXTURES}/sampleA.emapper.annotations" "${dir}/"
    sed 's/sampleA/sampleB/g' "${FIXTURES}/sampleA.emapper.annotations" > "${dir}/sampleB.emapper.annotations"
    cp "${GFF_FIXTURE}" "${dir}/sampleA.gff.gz"
    cp "${GFF_FIXTURE}" "${dir}/sampleB.gff.gz"
    cp "${FIXTURES}/sampleA.ags.tsv" "${dir}/"
    cp "${FIXTURES}/eggnog_versions.yml" "${dir}/"
}

# sampleB has no ags table on purpose: it exercises the TPM-only fallback
ASSEMBLY_ARGS=(--assembly-mode assembly
    --counts sampleA.featureCounts.txt sampleB.featureCounts.txt
    --annotations sampleA.emapper.annotations sampleB.emapper.annotations
    --gffs sampleA.gff.gz sampleB.gff.gz
    --ags sampleA.ags.tsv
    --eggnog-versions-yml eggnog_versions.yml --output-dir .)

#
# 1. Assembly-mode happy path
#
D="${WORK}/assembly"
seed_assembly_inputs "${D}"
NEXT_WORKDIR="${D}" expect_pass "assembly mode aggregates two samples" "${ASSEMBLY_ARGS[@]}"

for f in gene_annotations.tsv gene_abundance.tsv function_abundance.tsv \
         annotated_fraction.tsv ags_and_ge.tsv \
         function_wide_ko_tpm.tsv function_wide_cog_tpm.tsv \
         function_wide_ec_tpm.tsv function_wide_pfam_tpm.tsv function_wide_cazy_tpm.tsv \
         function_wide_cazy_dbcan_tpm.tsv \
         function_wide_ko_cpge.tsv function_wide_cog_cpge.tsv \
         function_wide_ec_cpge.tsv function_wide_pfam_cpge.tsv function_wide_cazy_cpge.tsv \
         function_wide_cazy_dbcan_cpge.tsv; do
    check "output ${f} exists" "[ -s '${D}/${f}' ]"
done

# No --dbcan in this run: the cazy_dbcan outputs are deterministic but empty
check "no --dbcan: cazy_dbcan matrices header-only, cazy_dbcan fraction 0" \
    "[ \$(wc -l < '${D}/function_wide_cazy_dbcan_tpm.tsv') -eq 1 ] && grep -qP 'sampleA\tcazy_dbcan\t2\t0\t0.0000\t0.0000' '${D}/annotated_fraction.tsv'"

# Acceptance 10.2: TPM sums to 1e6 per sample (both samples have reads)
check "TPM sums to 1e6 for both samples" \
    "awk -F'\t' 'NR>1 {sum[\$1]+=\$6} END {ok=1; for (s in sum) if (sum[s] < 999999.9 || sum[s] > 1000000.1) ok=0; exit !ok}' '${D}/gene_abundance.tsv'"

# Acceptance 10.2: every predicted gene appears, annotated or not (2 genes x 2 samples)
check "all 4 gene rows present incl. the count-0 gene" \
    "[ \$(tail -n +2 '${D}/gene_abundance.tsv' | wc -l) -eq 4 ] && grep -qP 'sampleB\tcontig_1_1\tassembly\t0\t' '${D}/gene_abundance.tsv'"

# Acceptance 10.2: multi-KO gene contributes full abundance to each KO
check "2-KO gene double-counts at full TPM (K00001 and K00002 both 750000)" \
    "grep -qP 'sampleA\tcontigs\teggnog-mapper\tko\tK00001\t\t750000.0000' '${D}/function_abundance.tsv' && grep -qP 'sampleA\tcontigs\teggnog-mapper\tko\tK00002\t\t750000.0000' '${D}/function_abundance.tsv'"

# Corrected 4.8 step 4: COG 'EG' splits per character
check "COG 'EG' yields separate E and G rows" \
    "grep -qP 'sampleA\tcontigs\teggnog-mapper\tcog\tE\t' '${D}/function_abundance.tsv' && grep -qP 'sampleA\tcontigs\teggnog-mapper\tcog\tG\t' '${D}/function_abundance.tsv' && ! grep -qP '\tcog\tEG\t' '${D}/function_abundance.tsv'"

# Q15 (2026-10-07): v3 PFAMs values are '<name>_<start>_<end>'; the
# coordinates are stripped and a repeated domain counts once per gene
# (fixture: MockPfam_10_95,MockPfam_150_290,Mock_dom_2_5_60)
check "pfam: coordinate suffix stripped, repeated domain counted once at full TPM" \
    "[ \$(grep -cP 'sampleA\tcontigs\teggnog-mapper\tpfam\tMockPfam\t\t750000.0000' '${D}/function_abundance.tsv') -eq 1 ] && [ \$(grep -cP '^sampleA\tcontig_1_1\t.*\teggnog_pfam\tMockPfam\t' '${D}/gene_annotations.tsv') -eq 1 ]"
check "pfam: name ending in digits kept whole (Mock_dom_2), no coordinate leftovers" \
    "grep -qP '\tpfam\tMock_dom_2\t' '${D}/function_abundance.tsv' && ! grep -qP '\tpfam\t[^\t]*_[0-9]+_[0-9]+\t' '${D}/function_abundance.tsv' && [ \$(wc -l < '${D}/function_wide_pfam_tpm.tsv') -eq 3 ]"

# '-' placeholder fields yield no rows (fixture EC is '-')
check "'-' EC field produces no ec rows" \
    "! grep -qP '\tec\t' '${D}/function_abundance.tsv' && [ \$(wc -l < '${D}/function_wide_ec_tpm.tsv') -eq 1 ]"

# Unannotated gene: absent from 5.1, present in 5.2
check "unannotated contig_1_2 absent from gene_annotations, present in gene_abundance" \
    "! grep -q 'contig_1_2' '${D}/gene_annotations.tsv' && grep -q 'contig_1_2' '${D}/gene_abundance.tsv'"

# Annotated fraction exact values (sampleA: 1/2 genes, 750000/1e6 TPM)
check "sampleA annotated fraction 0.5000 by count, 0.7500 by abundance" \
    "grep -qP 'sampleA\tany\t2\t1\t0.5000\t0.7500' '${D}/annotated_fraction.tsv'"

# Provenance columns carry the pinned versions
check "tool and DB versions in gene_annotations" \
    "grep -qP '\teggnog-mapper\t3.0.0-beta6\t7.0.0\$' '${D}/gene_annotations.tsv'"

# GFF ground truth: partial flag carried through (contig_1_1 partial=10)
check "partial flag carried from the GFF" \
    "grep -qP 'contig_1_1\tcontig_1\t1\t300\t\\+\t300\t10\t' '${D}/gene_annotations.tsv'"

# --- CPGE (T5, Section 6.2): hand-checked against the fixture GE of 2.0 ---
# RPK = count/(300/1000): 30 -> 100, 10 -> 33.333; cpge = RPK/2.0
check "sampleA cpge exact (50.000000 and 16.666667)" \
    "grep -qP 'sampleA\tcontig_1_1\tassembly\t30\t300\t750000.0000\t50.000000\$' '${D}/gene_abundance.tsv' && grep -qP 'sampleA\tcontig_1_2\tassembly\t10\t300\t250000.0000\t16.666667\$' '${D}/gene_abundance.tsv'"

# No ags table for sampleB -> empty cpge fields (TPM-only fallback), and TPM untouched
check "sampleB cpge empty, tpm unchanged" \
    "grep -qP 'sampleB\tcontig_1_1\tassembly\t0\t300\t0.0000\t\$' '${D}/gene_abundance.tsv' && grep -qP 'sampleB\tcontig_1_2\tassembly\t50\t300\t1000000.0000\t\$' '${D}/gene_abundance.tsv'"

# 5.3: the 2-KO gene double-counts its full CPGE into each KO
check "abundance_cpge 50.000000 for both sampleA KOs, empty for sampleB rows" \
    "grep -qP 'sampleA\tcontigs\teggnog-mapper\tko\tK00001\t\t750000.0000\t50.000000\t' '${D}/function_abundance.tsv' && grep -qP 'sampleA\tcontigs\teggnog-mapper\tko\tK00002\t\t750000.0000\t50.000000\t' '${D}/function_abundance.tsv' && grep -qP 'sampleB\tcontigs\teggnog-mapper\tko\tK00001\t\t0.0000\t\t' '${D}/function_abundance.tsv'"

# ags_and_ge.tsv: ok row with the fixture values, 'unavailable' warning row
check "ags_and_ge.tsv has sampleA ok and sampleB unavailable rows" \
    "grep -qP 'sampleA\t3000000\t2.0\t6000000\tok\$' '${D}/ags_and_ge.tsv' && grep -qP 'sampleB\t\t\t\tunavailable\$' '${D}/ags_and_ge.tsv'"

# Wide CPGE matrix: value column for sampleA, blank (not 0) column for sampleB
check "wide cpge matrix: sampleA values, sampleB blank" \
    "grep -qP '^K00001\t50.000000\t\$' '${D}/function_wide_ko_cpge.tsv' && head -1 '${D}/function_wide_ko_cpge.tsv' | grep -qP 'accession\tsampleA\tsampleB'"

# Backwards compatibility: no --ags at all still succeeds, all cpge empty
D="${WORK}/assembly_noags"
seed_assembly_inputs "${D}"
NEXT_WORKDIR="${D}" expect_pass "running without --ags succeeds (all-samples fallback)" \
    --assembly-mode assembly \
    --counts sampleA.featureCounts.txt sampleB.featureCounts.txt \
    --annotations sampleA.emapper.annotations sampleB.emapper.annotations \
    --gffs sampleA.gff.gz sampleB.gff.gz \
    --eggnog-versions-yml eggnog_versions.yml --output-dir .
check "no --ags: every ags_and_ge row unavailable, all cpge empty" \
    "[ \$(tail -n +2 '${D}/ags_and_ge.tsv' | grep -cP '\tunavailable\$') -eq 2 ] && [ -z \"\$(cut -f7 '${D}/gene_abundance.tsv' | tail -n +2 | tr -d '[:space:]')\" ]"

#
# 2. Coassembly mode (incl. zero-total sample)
#
D="${WORK}/coassembly"
mkdir -p "${D}"
cp "${FIXTURES}/sampleA.featureCounts.txt" "${FIXTURES}/sampleB.featureCounts.txt" \
   "${FIXTURES}/sampleZ.featureCounts.txt" "${D}/"
sed 's/sampleA/coassembly/g' "${FIXTURES}/sampleA.emapper.annotations" > "${D}/coassembly.emapper.annotations"
cp "${GFF_FIXTURE}" "${D}/coassembly.gff.gz"
cp "${FIXTURES}/sampleA.ags.tsv" "${D}/"
cp "${FIXTURES}/eggnog_versions.yml" "${D}/"
NEXT_WORKDIR="${D}" expect_pass "coassembly mode: shared gene set, per-sample counts" \
    --assembly-mode coassembly \
    --counts sampleA.featureCounts.txt sampleB.featureCounts.txt sampleZ.featureCounts.txt \
    --annotations coassembly.emapper.annotations \
    --gffs coassembly.gff.gz \
    --ags sampleA.ags.tsv \
    --eggnog-versions-yml eggnog_versions.yml --output-dir .

check "5.1 sample_id is 'coassembly' only" \
    "[ \"\$(tail -n +2 '${D}/gene_annotations.tsv' | cut -f1 | sort -u)\" = 'coassembly' ]"
check "wide matrix has all three sample columns" \
    "head -1 '${D}/function_wide_ko_tpm.tsv' | grep -qP 'accession\tsampleA\tsampleB\tsampleZ'"
check "zero-total sampleZ: tpm 0 rows, empty fraction_by_abundance" \
    "grep -qP 'sampleZ\tcontig_1_1\tcoassembly\t0\t300\t0.0000\t' '${D}/gene_abundance.tsv' && grep -qP 'sampleZ\tany\t2\t1\t0.5000\t\$' '${D}/annotated_fraction.tsv'"
check "coassembly: per-sample cpge (sampleA 50.000000, sampleB/Z empty)" \
    "grep -qP 'sampleA\tcontig_1_1\tcoassembly\t30\t300\t750000.0000\t50.000000\$' '${D}/gene_abundance.tsv' && grep -qP 'sampleB\tcontig_1_2\tcoassembly\t50\t300\t1000000.0000\t\$' '${D}/gene_abundance.tsv' && grep -qP 'sampleZ\t\t\t\tunavailable\$' '${D}/ags_and_ge.tsv'"

#
# 3. Empty gene set (empty-contig sample shape)
#
D="${WORK}/empty"
mkdir -p "${D}"
printf '# Program:featureCounts v2.1.1; empty gene set\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tsampleE_all_reads.bam\n' \
    > "${D}/sampleE.featureCounts.txt"
# Header-only annotations (an emapper run over zero proteins)
grep '^#' "${FIXTURES}/sampleA.emapper.annotations" > "${D}/sampleE.emapper.annotations"
cp "${EMPTY_GFF_FIXTURE}" "${D}/sampleE.gff.gz"
cp "${FIXTURES}/eggnog_versions.yml" "${D}/"
NEXT_WORKDIR="${D}" expect_pass "empty gene set aggregates without crashing" \
    --assembly-mode assembly \
    --counts sampleE.featureCounts.txt --annotations sampleE.emapper.annotations \
    --gffs sampleE.gff.gz --eggnog-versions-yml eggnog_versions.yml --output-dir .

check "empty sample: zero abundance rows, summary row with empty fractions" \
    "[ \$(tail -n +2 '${D}/gene_abundance.tsv' | wc -l) -eq 0 ] && grep -qP 'sampleE\tany\t0\t0\t\t\$' '${D}/annotated_fraction.tsv' && grep -qP 'sampleE\t\t\t\tunavailable\$' '${D}/ags_and_ge.tsv'"

#
# 4. Version-aware layout guard: deliberate loud failures (acceptance 10.2)
#
seed_negative() {
    # seed_negative <dir> — single-sample assembly inputs to mutate
    local dir="$1"
    mkdir -p "${dir}"
    cp "${FIXTURES}/sampleA.featureCounts.txt" "${dir}/"
    cp "${FIXTURES}/sampleA.emapper.annotations" "${dir}/"
    cp "${GFF_FIXTURE}" "${dir}/sampleA.gff.gz"
    cp "${FIXTURES}/eggnog_versions.yml" "${dir}/"
}
NEG_ARGS=(--assembly-mode assembly --counts sampleA.featureCounts.txt
    --annotations sampleA.emapper.annotations --gffs sampleA.gff.gz
    --eggnog-versions-yml eggnog_versions.yml --output-dir .)

D="${WORK}/neg_renamed"; seed_negative "${D}"
sed -i 's/\tKEGG_ko\t/\tKEGG_KO_RENAMED\t/' "${D}/sampleA.emapper.annotations"
NEXT_WORKDIR="${D}" expect_fail "renamed column in annotations header fails" "${NEG_ARGS[@]}"
check_grep "  ...with a layout-drift message" "../neg_renamed.log" "does not match the verified"

D="${WORK}/neg_21col"; seed_negative "${D}"
sed -i 's/\tannotation_confidence$//; s/\th--h------h-h$//' "${D}/sampleA.emapper.annotations"
NEXT_WORKDIR="${D}" expect_fail "21-column annotations header fails" "${NEG_ARGS[@]}"

D="${WORK}/neg_23col"; seed_negative "${D}"
sed -i 's/\tannotation_confidence$/\tannotation_confidence\tmd5/' "${D}/sampleA.emapper.annotations"
NEXT_WORKDIR="${D}" expect_fail "23-column (--md5) annotations header fails" "${NEG_ARGS[@]}"

D="${WORK}/neg_unknown_version"; seed_negative "${D}"
sed -i 's/eggnog-mapper: 3.0.0-beta6/eggnog-mapper: 9.9.9/' "${D}/eggnog_versions.yml"
NEXT_WORKDIR="${D}" expect_fail "unverified emapper version fails" "${NEG_ARGS[@]}"
check_grep "  ...naming the known-layout set" "../neg_unknown_version.log" "not among the layouts"

D="${WORK}/neg_missing_version"; seed_negative "${D}"
grep -v 'eggnog-mapper:' "${FIXTURES}/eggnog_versions.yml" > "${D}/eggnog_versions.yml"
NEXT_WORKDIR="${D}" expect_fail "missing emapper version in versions.yml fails" "${NEG_ARGS[@]}"

D="${WORK}/neg_pfam_no_coords"; seed_negative "${D}"
sed -i 's/MockPfam_10_95,MockPfam_150_290,Mock_dom_2_5_60/MockPfam/' "${D}/sampleA.emapper.annotations"
NEXT_WORKDIR="${D}" expect_fail "PFAMs value without the v3 coordinate suffix fails (Q15)" "${NEG_ARGS[@]}"
check_grep "  ...naming the offending PFAMs value" "../neg_pfam_no_coords.log" "PFAMs value .MockPfam. is not"

D="${WORK}/neg_header_disagrees"; seed_negative "${D}"
sed -i 's/^## emapper-3.0.0-beta6$/## emapper-2.1.15/' "${D}/sampleA.emapper.annotations"
NEXT_WORKDIR="${D}" expect_fail "'## emapper-' line disagreeing with versions.yml fails" "${NEG_ARGS[@]}"

D="${WORK}/neg_id_mismatch"; seed_negative "${D}"
printf 'contig_1_99\tcontig_1\t950\t999\t+\t50\t7\n' >> "${D}/sampleA.featureCounts.txt"
NEXT_WORKDIR="${D}" expect_fail "counts/GFF gene id mismatch fails" "${NEG_ARGS[@]}"

D="${WORK}/neg_length_mismatch"; seed_negative "${D}"
sed -i 's/contig_1_1\tcontig_1\t1\t300\t+\t300\t30/contig_1_1\tcontig_1\t1\t300\t+\t299\t30/' "${D}/sampleA.featureCounts.txt"
NEXT_WORKDIR="${D}" expect_fail "featureCounts Length vs GFF length mismatch fails" "${NEG_ARGS[@]}"

#
# 5. AGS guard: deliberate loud failures (T5)
#
D="${WORK}/neg_ags_header"; seed_negative "${D}"
printf 'sample\tags\tge\ttb\nsampleA\t3000000\t2.0\t6000000\n' > "${D}/sampleA.ags.tsv"
NEXT_WORKDIR="${D}" expect_fail "wrong ags.tsv header fails" "${NEG_ARGS[@]}" --ags sampleA.ags.tsv
check_grep "  ...naming the expected header" "../neg_ags_header.log" "unexpected ags.tsv header"

D="${WORK}/neg_ags_unknown"; seed_negative "${D}"
sed 's/^sampleA\t/sampleQ\t/' "${FIXTURES}/sampleA.ags.tsv" > "${D}/sampleQ.ags.tsv"
NEXT_WORKDIR="${D}" expect_fail "ags table for a sample not in --counts fails" "${NEG_ARGS[@]}" --ags sampleQ.ags.tsv

D="${WORK}/neg_ags_zero_ge"; seed_negative "${D}"
sed 's/\t2.0\t/\t0\t/' "${FIXTURES}/sampleA.ags.tsv" > "${D}/sampleA.ags.tsv"
NEXT_WORKDIR="${D}" expect_fail "non-positive genome_equivalents fails" "${NEG_ARGS[@]}" --ags sampleA.ags.tsv
check_grep "  ...naming the offending field" "../neg_ags_zero_ge.log" "genome_equivalents must be a positive number"

D="${WORK}/neg_ags_id_mismatch"; seed_negative "${D}"
cp "${FIXTURES}/sampleA.ags.tsv" "${D}/sampleM.ags.tsv"
NEXT_WORKDIR="${D}" expect_fail "ags filename/embedded sample_id mismatch fails" "${NEG_ARGS[@]}" --ags sampleM.ags.tsv

#
# 6. dbCAN CAZy branch (T6): happy paths
#
seed_dbcan_inputs() {
    # seed_assembly_inputs plus per-sample dbCAN overviews and versions.yml
    local dir="$1"
    seed_assembly_inputs "${dir}"
    cp "${FIXTURES}/sampleA.overview.tsv" "${dir}/"
    cp "${FIXTURES}/sampleA.overview.tsv" "${dir}/sampleB.overview.tsv"
    cp "${FIXTURES}/dbcan_versions.yml" "${dir}/"
}
DBCAN_ARGS=("${ASSEMBLY_ARGS[@]}"
    --dbcan sampleA.overview.tsv sampleB.overview.tsv
    --dbcan-versions-yml dbcan_versions.yml)

D="${WORK}/dbcan_assembly"
seed_dbcan_inputs "${D}"
NEXT_WORKDIR="${D}" expect_pass "assembly mode with dbCAN overviews" "${DBCAN_ARGS[@]}"

# 5.1: db=dbcan rows with per-row run_dbcan provenance; empty evalue/score
check "dbcan 5.1 row with run_dbcan provenance (GH5_4 on contig_1_1)" \
    "grep -qP 'sampleA\tcontig_1_1\tcontig_1\t1\t300\t\\+\t300\t10\tdbcan\tGH5_4\t\t\t\trun_dbcan\t5.2.9\tdb_v5-2-9_5-5-2026\$' '${D}/gene_annotations.tsv'"
check "eggNOG rows keep their own provenance next to dbcan rows" \
    "grep -qP '\teggnog_cazy\tGT2\t.*\teggnog-mapper\t3.0.0-beta6\t7.0.0\$' '${D}/gene_annotations.tsv'"

# Consensus 'recommended' (default): the tool's >=2-tools column {GH5_4, CBM6};
# the hmm-only GH5 call must NOT appear
check "recommended consensus: GH5_4 and CBM6 in, hmm-only GH5 out" \
    "grep -qP '\tdbcan\tCBM6\t' '${D}/gene_annotations.tsv' && ! grep -qP '\tdbcan\tGH5\t' '${D}/gene_annotations.tsv'"

# 5.3: backend=run_dbcan rows mirror the gene's TPM/CPGE; eggNOG cazy intact
check "dbcan 5.3 row carries the gene's TPM and CPGE" \
    "grep -qP 'sampleA\tcontigs\trun_dbcan\tcazy\tGH5_4\t\t750000.0000\t50.000000\t' '${D}/function_abundance.tsv' && grep -qP 'sampleA\tcontigs\teggnog-mapper\tcazy\tGT2\t' '${D}/function_abundance.tsv'"

# Owner decision (4.3): separate wide matrices, never merged
check "cazy and cazy_dbcan wide matrices stay separate" \
    "grep -qP '^GH5_4\t' '${D}/function_wide_cazy_dbcan_tpm.tsv' && ! grep -q 'GT2' '${D}/function_wide_cazy_dbcan_tpm.tsv' && grep -qP '^GT2\t' '${D}/function_wide_cazy_tpm.tsv' && ! grep -q 'GH5_4' '${D}/function_wide_cazy_tpm.tsv'"
check "cazy_dbcan cpge matrix: sampleA value, AGS-less sampleB blank" \
    "grep -qP '^GH5_4\t50.000000\t\$' '${D}/function_wide_cazy_dbcan_cpge.tsv'"

# Fractions: 'cazy' counts either backend, 'cazy_dbcan' counts run_dbcan alone
check "cazy and cazy_dbcan fractions for sampleA (1/2 genes, 0.75 by TPM)" \
    "grep -qP 'sampleA\tcazy\t2\t1\t0.5000\t0.7500' '${D}/annotated_fraction.tsv' && grep -qP 'sampleA\tcazy_dbcan\t2\t1\t0.5000\t0.7500' '${D}/annotated_fraction.tsv'"

# Consensus 'any': union of the per-tool columns, so GH5 appears
D="${WORK}/dbcan_any"
seed_dbcan_inputs "${D}"
NEXT_WORKDIR="${D}" expect_pass "'any' consensus takes the per-tool union" \
    "${DBCAN_ARGS[@]}" --dbcan-consensus any
check "any: range-stripped hmm-only GH5 row present" \
    "grep -qP '\tdbcan\tGH5\t' '${D}/gene_annotations.tsv' && grep -qP '\tdbcan\tGH5_4\t' '${D}/gene_annotations.tsv'"

# Coassembly: one shared overview, per-sample 5.3 rows
D="${WORK}/dbcan_coassembly"
mkdir -p "${D}"
cp "${FIXTURES}/sampleA.featureCounts.txt" "${FIXTURES}/sampleB.featureCounts.txt" "${D}/"
cp "${FIXTURES}/coassembly.emapper.annotations" "${D}/"
cp "${GFF_FIXTURE}" "${D}/coassembly.gff.gz"
cp "${FIXTURES}/coassembly.overview.tsv" "${D}/"
cp "${FIXTURES}/eggnog_versions.yml" "${FIXTURES}/dbcan_versions.yml" "${D}/"
NEXT_WORKDIR="${D}" expect_pass "coassembly mode with one shared dbCAN overview" \
    --assembly-mode coassembly \
    --counts sampleA.featureCounts.txt sampleB.featureCounts.txt \
    --annotations coassembly.emapper.annotations \
    --gffs coassembly.gff.gz \
    --dbcan coassembly.overview.tsv \
    --dbcan-versions-yml dbcan_versions.yml \
    --eggnog-versions-yml eggnog_versions.yml --output-dir .
check "coassembly: dbcan 5.1 keyed 'coassembly', per-sample 5.3 rows" \
    "grep -qP 'coassembly\tcontig_1_1\tcontig_1\t1\t300\t\\+\t300\t10\tdbcan\tGH5_4\t' '${D}/gene_annotations.tsv' && grep -qP 'sampleA\tcontigs\trun_dbcan\tcazy\tGH5_4\t' '${D}/function_abundance.tsv' && grep -qP 'sampleB\tcontigs\trun_dbcan\tcazy\tGH5_4\t' '${D}/function_abundance.tsv'"

#
# 7. dbCAN guard: deliberate loud failures (acceptance 10.2)
#
seed_dbcan_negative() {
    # seed_negative plus a single-sample overview and versions.yml to mutate
    local dir="$1"
    seed_negative "${dir}"
    cp "${FIXTURES}/sampleA.overview.tsv" "${dir}/"
    cp "${FIXTURES}/dbcan_versions.yml" "${dir}/"
}
DBCAN_NEG_ARGS=("${NEG_ARGS[@]}" --dbcan sampleA.overview.tsv
    --dbcan-versions-yml dbcan_versions.yml)

D="${WORK}/neg_dbcan_header"; seed_dbcan_negative "${D}"
sed -i 's/^Gene ID\t/GeneID\t/' "${D}/sampleA.overview.tsv"
NEXT_WORKDIR="${D}" expect_fail "corrupt overview header fails" "${DBCAN_NEG_ARGS[@]}"
check_grep "  ...with a layout-drift message" "../neg_dbcan_header.log" "does not match the verified"

D="${WORK}/neg_dbcan_row"; seed_dbcan_negative "${D}"
printf 'contig_1_1\tGH1\tbad\n' >> "${D}/sampleA.overview.tsv"
NEXT_WORKDIR="${D}" expect_fail "overview row with the wrong field count fails" "${DBCAN_NEG_ARGS[@]}"

D="${WORK}/neg_dbcan_version"; seed_dbcan_negative "${D}"
sed -i 's/run_dbcan: 5.2.9/run_dbcan: 9.9.9/' "${D}/dbcan_versions.yml"
NEXT_WORKDIR="${D}" expect_fail "unverified run_dbcan version fails" "${DBCAN_NEG_ARGS[@]}"
check_grep "  ...naming the known-layout set" "../neg_dbcan_version.log" "not among the layouts"

D="${WORK}/neg_dbcan_noversions"; seed_dbcan_negative "${D}"
NEXT_WORKDIR="${D}" expect_fail "--dbcan without --dbcan-versions-yml fails" \
    "${NEG_ARGS[@]}" --dbcan sampleA.overview.tsv
check_grep "  ...naming the missing input" "../neg_dbcan_noversions.log" "without --dbcan-versions-yml"

D="${WORK}/neg_dbcan_stray"; seed_dbcan_negative "${D}"
printf 'contig_1_99\t-\tGH1(1-50)\t-\t-\t1\tGH1\t-\n' >> "${D}/sampleA.overview.tsv"
NEXT_WORKDIR="${D}" expect_fail "dbCAN gene id not in the GFF fails" "${DBCAN_NEG_ARGS[@]}"
check_grep "  ...naming the stray id" "../neg_dbcan_stray.log" "dbCAN gene ids not present"

D="${WORK}/neg_dbcan_sample"; seed_dbcan_negative "${D}"
mv "${D}/sampleA.overview.tsv" "${D}/sampleQ.overview.tsv"
NEXT_WORKDIR="${D}" expect_fail "dbCAN sample set not matching --counts fails" \
    "${NEG_ARGS[@]}" --dbcan sampleQ.overview.tsv --dbcan-versions-yml dbcan_versions.yml

D="${WORK}/neg_dbcan_duplicate"; seed_dbcan_negative "${D}"
NEXT_WORKDIR="${D}" expect_fail "duplicate overview for one sample fails" \
    "${NEG_ARGS[@]}" --dbcan sampleA.overview.tsv sampleA.overview.tsv \
    --dbcan-versions-yml dbcan_versions.yml

#
# Summary
#
echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
