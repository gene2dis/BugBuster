#!/bin/bash
#
# Direct tests for bin/aggregate_functions.py (functional annotation T4,
# design doc Sections 4.8, 5, 6.1, 10.2). Asserts on the committed fixtures
# in tests/data/functional/ (regenerate with make_functional_fixtures.sh):
#   - assembly mode: TPM sums to 1e6 per sample; every predicted gene appears
#     in gene_abundance including count-0 genes; a gene with two KOs
#     contributes its full abundance to each (intentional double-counting);
#     multi-letter COG splits per character; '-' fields yield no rows;
#     unannotated genes are absent from 5.1 but present in 5.2; annotated
#     fractions by count and abundance are exact
#   - coassembly mode: one shared annotations/GFF, per-sample counts;
#     5.1 sample_id is 'coassembly'; zero-total sample gets tpm 0 and an
#     empty fraction_by_abundance
#   - empty gene set: header-only inputs aggregate without crashing
#   - the version-aware layout guard fails loudly (10.2 mandates testing
#     this deliberately): renamed / 21-column / 23-column headers, unknown
#     or missing emapper version, a disagreeing '## emapper-' line, and
#     mis-paired gene id sets all exit non-zero
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
    cp "${FIXTURES}/eggnog_versions.yml" "${dir}/"
}

ASSEMBLY_ARGS=(--assembly-mode assembly
    --counts sampleA.featureCounts.txt sampleB.featureCounts.txt
    --annotations sampleA.emapper.annotations sampleB.emapper.annotations
    --gffs sampleA.gff.gz sampleB.gff.gz
    --eggnog-versions-yml eggnog_versions.yml --output-dir .)

#
# 1. Assembly-mode happy path
#
D="${WORK}/assembly"
seed_assembly_inputs "${D}"
NEXT_WORKDIR="${D}" expect_pass "assembly mode aggregates two samples" "${ASSEMBLY_ARGS[@]}"

for f in gene_annotations.tsv gene_abundance.tsv function_abundance.tsv \
         annotated_fraction.tsv function_wide_ko_tpm.tsv function_wide_cog_tpm.tsv \
         function_wide_ec_tpm.tsv function_wide_pfam_tpm.tsv function_wide_cazy_tpm.tsv; do
    check "output ${f} exists" "[ -s '${D}/${f}' ]"
done

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

#
# 2. Coassembly mode (incl. zero-total sample)
#
D="${WORK}/coassembly"
mkdir -p "${D}"
cp "${FIXTURES}/sampleA.featureCounts.txt" "${FIXTURES}/sampleB.featureCounts.txt" \
   "${FIXTURES}/sampleZ.featureCounts.txt" "${D}/"
sed 's/sampleA/coassembly/g' "${FIXTURES}/sampleA.emapper.annotations" > "${D}/coassembly.emapper.annotations"
cp "${GFF_FIXTURE}" "${D}/coassembly.gff.gz"
cp "${FIXTURES}/eggnog_versions.yml" "${D}/"
NEXT_WORKDIR="${D}" expect_pass "coassembly mode: shared gene set, per-sample counts" \
    --assembly-mode coassembly \
    --counts sampleA.featureCounts.txt sampleB.featureCounts.txt sampleZ.featureCounts.txt \
    --annotations coassembly.emapper.annotations \
    --gffs coassembly.gff.gz \
    --eggnog-versions-yml eggnog_versions.yml --output-dir .

check "5.1 sample_id is 'coassembly' only" \
    "[ \"\$(tail -n +2 '${D}/gene_annotations.tsv' | cut -f1 | sort -u)\" = 'coassembly' ]"
check "wide matrix has all three sample columns" \
    "head -1 '${D}/function_wide_ko_tpm.tsv' | grep -qP 'accession\tsampleA\tsampleB\tsampleZ'"
check "zero-total sampleZ: tpm 0 rows, empty fraction_by_abundance" \
    "grep -qP 'sampleZ\tcontig_1_1\tcoassembly\t0\t300\t0.0000\t' '${D}/gene_abundance.tsv' && grep -qP 'sampleZ\tany\t2\t1\t0.5000\t\$' '${D}/annotated_fraction.tsv'"

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
    "[ \$(tail -n +2 '${D}/gene_abundance.tsv' | wc -l) -eq 0 ] && grep -qP 'sampleE\tany\t0\t0\t\t\$' '${D}/annotated_fraction.tsv'"

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
# Summary
#
echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
