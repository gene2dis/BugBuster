#!/bin/bash
#
# Direct tests for the read-level functional branch scripts (design doc
# Sections 4.6, 4.6.2, 4.6.3, 5.3, 6.2, 10.2; tasks T8a - Woltka backend and
# T8b - the SUPER-FOCUS backend's aggregation), each run
# inside the exact image its module pins:
#
#   bin/woltka_function_profile.py   (woltka image, WOLTKA_CLASSIFY)
#     - real `woltka classify` on the committed fixture alignment
#       (tests/data/woltka/wtest.wol.sam.gz, made by make_woltka_fixtures.sh)
#       with the module's fixed flags: mates count separately (4 pairs -> 8),
#       ambiguous reads divide 1/k across the two shared-region ORFs
#     - term-set composition is de-duplicated per ORF (the reason woltka
#       collapse is not used): a 3x-repeated Pfam domain counts once, two
#       KOs mapping to one EC count the ORF once, a MetaCyc pathway reached
#       through several enzrxn/reaction paths counts the ORF once
#     - RPK = count / (length / 1000), names attached, annotated fractions
#     - Q14 vocabularies: COG ids -> category letters via the COG table
#       (tests/data/woltka/cog.def.tab), de-duplicated per ORF AFTER mapping
#       (an ORF reaching E through two COGs counts once); the skip-listed
#       WoLr2 defect COG:1140 is skipped with a warning; pfam rows carry the
#       Pfam NAME with the versioned accession as description
#     - --uniq: ambiguous mates move to the unassigned table
#     - header-only (empty-alignment) profile -> header-only functions,
#       zero summary rows with empty fractions
#     - guards fail loudly: a demultiplexed (3-column) profile, a profile ORF
#       missing from length.map, a missing database file, a COG id neither
#       in the table nor skip-listed, no --cog-def, a malformed COG table, a
#       Pfam accession without a name, a Pfam name shared by two accessions
#   bin/aggregate_read_functions.py  (pandas image, AGGREGATE_READ_FUNCTIONS)
#     - Section 5.3 schema with source=reads, backend=woltka, empty
#       abundance_tpm, abundance_native = read counts, native_unit=reads
#     - CPGE = RPK / genome_equivalents exactly; blank (never 0) for samples
#       without AGS, with cpge_status 'unavailable' in read_sample_summary
#     - wide native/cpge matrices for all five ontologies, always written
#     - version-aware guard (unknown or missing woltka version), header,
#       sample-set, ontology, duplicate and ags guards all fail loudly
#     - superfocus backend (committed tests/data/superfocus composer
#       outputs): seed_level1..3 only, CPGE blank / not_applicable with AGS
#       still reported, only the three native wide matrices, and its own
#       version / rpk / ontology / --unassigned guards
#       (the SUPER-FOCUS composer itself: tests/bin/test_superfocus_scripts.sh)
#
# Requirements: docker (present on GitHub ubuntu-latest runners).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
FIX="${REPO_DIR}/tests/data/woltka"

# Exact images pinned in the modules, so the test cannot drift from them
WOLTKA_IMG=$(grep -o "'quay.io/biocontainers/woltka:[^']*'" "${REPO_DIR}/modules/local/woltka_classify/main.nf" \
    | tr -d "'" | head -1)
REPORT_IMG=$(grep -o "'[^']*community.wave.seqera.io/library/python_pandas[^']*'" "${REPO_DIR}/modules/local/aggregate_read_functions/main.nf" \
    | tr -d "'" | grep -v '^oras://' | head -1)
if [ -z "${WOLTKA_IMG}" ] || [ -z "${REPORT_IMG}" ]; then
    echo "ERROR: could not extract pinned container images from the modules"
    exit 1
fi

echo "=== read-level functional script test suite ==="
echo "Woltka image: ${WOLTKA_IMG}"
echo "Report image: ${REPORT_IMG}"
echo ""

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0

in_woltka() {
    # in_woltka <workdir> <shell command> - fixture mounted read-only at /fix
    local workdir="$1"; shift
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" -v "${FIX}:/fix:ro" \
        -v "${workdir}:/w" -w /w "${WOLTKA_IMG}" sh -c "$*"
}

in_report() {
    local workdir="$1"; shift
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e MPLCONFIGDIR=/tmp \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" \
        -v "${workdir}:/w" -w /w "${REPORT_IMG}" \
        python3 /pipeline_bin/aggregate_read_functions.py "$@"
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

expect_ok() {
    local desc="$1" log="$2"; shift 2
    if "$@" > "${log}" 2>&1; then
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${desc} - expected success:"
        tail -5 "${log}" | sed 's/^/    /'
        FAIL=$((FAIL + 1))
    fi
}

expect_err() {
    # expect_err <desc> <log> <expected-message-fragment> <command...>
    local desc="$1" log="$2" msg="$3"; shift 3
    if "$@" > "${log}" 2>&1; then
        echo "✗ ${desc} - expected failure, got success"
        FAIL=$((FAIL + 1))
    elif grep -q -- "${msg}" "${log}"; then
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${desc} - failed, but without '${msg}':"
        tail -3 "${log}" | sed 's/^/    /'
        FAIL=$((FAIL + 1))
    fi
}

# value <tsv> <awk-condition> <column-number>
value() { awk -F'\t' -v c="$3" "$2 { print \$c }" "$1"; }

CLASSIFY='woltka classify --input /fix/wtest.wol.sam.gz --coords /fix/db/proteins/coords.txt.xz --no-demux --digits 6 --unassigned --to-tsv'
PROFILE='python3 /pipeline_bin/woltka_function_profile.py --db /fix/db --cog-def /fix/cog.def.tab'

#
# woltka_function_profile.py
#
echo "--- woltka_function_profile.py (woltka image) ---"
D="${TMP_DIR}/profile"; mkdir -p "${D}"
expect_ok "classify + compose on the fixture alignment" "${D}/run.log" \
    in_woltka "${D}" "${CLASSIFY} --output orf.tsv && ${PROFILE} --profile orf.tsv --sample-id wtest --prefix wtest"
F="${D}/wtest.woltka_functions.tsv"
SUM="${D}/wtest.woltka_summary.tsv"
check "mates count separately: 4 pairs in G000000001_2 -> 8 reads" \
    "[ \"\$(value ${D}/orf.tsv '\$1==\"G000000001_2\"' 2)\" = '8.0' ]"
check "ambiguous shared-region pairs divide 1/k: 3 per shared ORF" \
    "[ \"\$(value ${D}/orf.tsv '\$1==\"G000000001_1\"' 2)\" = '3.0' ] && [ \"\$(value ${D}/orf.tsv '\$1==\"G000000002_1\"' 2)\" = '3.0' ]"
check "header matches the composer layout" \
    "[ \"\$(head -1 ${F})\" = \"\$(printf 'sample_id\tontology\taccession\tdescription\tcount\trpk')\" ]"
check "K00001 = 6 reads (two ORFs), named from ko_name.txt" \
    "[ \"\$(value ${F} '\$2==\"ko\" && \$3==\"K00001\"' 5)\" = '6' ] && grep -q 'alcohol dehydrogenase' ${F}"
check "EC 1.1.1.1 counts the 2-KO ORF once (6, not woltka-collapse's 9)" \
    "[ \"\$(value ${F} '\$2==\"ec\" && \$3==\"1.1.1.1\"' 5)\" = '6' ]"
check "repeated Pfam domain counts once (7tm_2 = 3, not 9)" \
    "[ \"\$(value ${F} '\$2==\"pfam\" && \$3==\"7tm_2\"' 5)\" = '3' ]"
check "pfam rows: Pfam name as accession, versioned accession as description" \
    "[ \"\$(value ${F} '\$2==\"pfam\"' 3 | tr '\n' ' ')\" = '7tm_1 7tm_2 7tm_3 ' ] && [ \"\$(value ${F} '\$2==\"pfam\" && \$3==\"7tm_1\"' 4)\" = 'PF00001.1' ]"
check "multi-path MetaCyc pathway counts each ORF once (PWY-1 = 7, not 10)" \
    "[ \"\$(value ${F} '\$2==\"metacyc\" && \$3==\"PWY-1\"' 5)\" = '7' ]"
check "cog rows are category letters only (C E H)" \
    "[ \"\$(value ${F} '\$2==\"cog\"' 3 | tr '\n' ' ')\" = 'C E H ' ]"
check "COG letters de-duplicated per ORF after mapping (E = 14, not 22)" \
    "[ \"\$(value ${F} '\$2==\"cog\" && \$3==\"E\"' 5)\" = '14' ] && [ \"\$(value ${F} '\$2==\"cog\" && \$3==\"H\"' 5)\" = '14' ]"
check "COG category whitespace stripped ('C ' -> C = 3)" \
    "[ \"\$(value ${F} '\$2==\"cog\" && \$3==\"C\"' 5)\" = '3' ]"
check "skip-listed WoLr2 COG defect: warned, never an accession" \
    "grep -q 'known WoLr2 defects: COG:1140' ${D}/run.log && ! grep -q 'COG' ${F}"
check "RPK = count / kb (K00003: 8 reads / 1.2 kb = 6.6667)" \
    "value ${F} '\$2==\"ko\" && \$3==\"K00003\"' 6 | grep -q '^6\.66666666'"
check "ec and cog carry no description" \
    "[ -z \"\$(value ${F} '\$2==\"ec\" || \$2==\"cog\"' 4 | tr -d '\n')\" ]"
check "matches the committed fixture output" \
    "cmp -s ${F} ${FIX}/wtest.woltka_functions.tsv && cmp -s ${SUM} ${FIX}/wtest.woltka_summary.tsv"
check "summary: any 18/20 annotated (G000000001_3 unannotated)" \
    "[ \"\$(value ${SUM} '\$2==\"any\"' 3)\" = '20' ] && [ \"\$(value ${SUM} '\$2==\"any\"' 5)\" = '0.900000' ]"
check "summary: metacyc 7/20" \
    "[ \"\$(value ${SUM} '\$2==\"metacyc\"' 4)\" = '7' ]"
check "default mode: no reads unassigned" \
    "[ \"\$(value ${D}/wtest.woltka_unassigned.tsv 'NR==2' 2)\" = '0' ]"

U="${TMP_DIR}/uniq"; mkdir -p "${U}"
expect_ok "classify --uniq + compose" "${U}/run.log" \
    in_woltka "${U}" "${CLASSIFY} --uniq --output orf.tsv && ${PROFILE} --profile orf.tsv --sample-id wtest --prefix wtest"
check "--uniq: the 6 ambiguous mates are unassigned" \
    "[ \"\$(value ${U}/wtest.woltka_unassigned.tsv 'NR==2' 2)\" = '6' ]"
check "--uniq: K00001 drops to 0 rows, K00003 unchanged" \
    "! grep -q 'K00001' ${U}/wtest.woltka_functions.tsv && [ \"\$(value ${U}/wtest.woltka_functions.tsv '\$3==\"K00003\"' 5)\" = '8' ]"

E="${TMP_DIR}/empty"; mkdir -p "${E}"
printf '#FeatureID\tempty\n' > "${E}/orf.tsv"
expect_ok "header-only (no alignments) profile" "${E}/run.log" \
    in_woltka "${E}" "${PROFILE} --profile orf.tsv --sample-id empty --prefix empty"
check "empty: header-only functions table" "[ \"\$(wc -l < ${E}/empty.woltka_functions.tsv)\" -eq 1 ]"
check "empty: six zero summary rows with empty fractions" \
    "[ \"\$(awk -F'\t' 'NR>1 && \$3==\"0\" && \$5==\"\"' ${E}/empty.woltka_summary.tsv | wc -l)\" -eq 6 ]"

N="${TMP_DIR}/neg"; mkdir -p "${N}"
printf '#FeatureID\ts\tx\nG000000001_1\t3.0\t0\n' > "${N}/demux.tsv"
printf '#FeatureID\tz\nG999999999_1\t2.0\n' > "${N}/unknown_orf.tsv"
cp "${D}/orf.tsv" "${N}/orf.tsv"
expect_err "guard: demultiplexed 3-column profile fails" "${N}/demux.log" "unexpected Woltka profile header" \
    in_woltka "${N}" "${PROFILE} --profile demux.tsv --sample-id s --prefix s"
expect_err "guard: profile ORF missing from length.map fails" "${N}/orf.log" "have no entry in proteins/length.map.xz" \
    in_woltka "${N}" "${PROFILE} --profile unknown_orf.tsv --sample-id z --prefix z"
expect_err "guard: missing database file fails" "${N}/db.log" "Woltka database file missing" \
    in_woltka "${N}" "mkdir -p partial && cp -r /fix/db/proteins partial/ && python3 /pipeline_bin/woltka_function_profile.py --db partial --cog-def /fix/cog.def.tab --profile orf.tsv --sample-id w --prefix w"
expect_err "guard: COG id neither in the table nor skip-listed fails" "${N}/cog.log" "is not a known WoLr2" \
    in_woltka "${N}" "rm -rf badcog && cp -r /fix/db badcog && printf 'K09999\tCOG9999\n' >> badcog/function/kegg/ko-to-cog.map && python3 /pipeline_bin/woltka_function_profile.py --db badcog --cog-def /fix/cog.def.tab --profile orf.tsv --sample-id w --prefix w"
expect_err "guard: --cog-def is required" "${N}/nocog.log" "the following arguments are required: --cog-def" \
    in_woltka "${N}" "python3 /pipeline_bin/woltka_function_profile.py --db /fix/db --profile orf.tsv --sample-id w --prefix w"
expect_err "guard: malformed COG table fails" "${N}/badtab.log" "malformed COG definition" \
    in_woltka "${N}" "printf 'COG0001\t1\tbad\n' > bad.def.tab && python3 /pipeline_bin/woltka_function_profile.py --db /fix/db --cog-def bad.def.tab --profile orf.tsv --sample-id w --prefix w"
expect_err "guard: Pfam accession without a name fails" "${N}/pfname.log" "has no entry in function/pfam/pfam_name.txt" \
    in_woltka "${N}" "rm -rf nopfam && cp -r /fix/db nopfam && grep -v PF00003.1 /fix/db/function/pfam/pfam_name.txt > nopfam/function/pfam/pfam_name.txt && python3 /pipeline_bin/woltka_function_profile.py --db nopfam --cog-def /fix/cog.def.tab --profile orf.tsv --sample-id w --prefix w"
expect_err "guard: Pfam name shared by two accessions fails" "${N}/pfdup.log" "is shared by" \
    in_woltka "${N}" "rm -rf duppfam && cp -r /fix/db duppfam && sed 's/7tm_3/7tm_2/' /fix/db/function/pfam/pfam_name.txt > duppfam/function/pfam/pfam_name.txt && python3 /pipeline_bin/woltka_function_profile.py --db duppfam --cog-def /fix/cog.def.tab --profile orf.tsv --sample-id w --prefix w"

#
# aggregate_read_functions.py
#
echo "--- aggregate_read_functions.py (report image) ---"
A="${TMP_DIR}/agg"; mkdir -p "${A}/out"
cp "${D}"/wtest.woltka_*.tsv "${A}/"
for f in "${U}"/wtest.woltka_*.tsv; do
    sed 's/^wtest\t/wuniq\t/' "${f}" > "${A}/wuniq.${f##*/wtest.}"
done
printf 'sample_id\taverage_genome_size_bp\tgenome_equivalents\ttotal_bases\nwtest\t3000000\t2.0\t6000000\n' > "${A}/wtest.ags.tsv"
printf '"WOLTKA_CLASSIFY":\n    woltka: 0.1.7\n    wol_db: WoLr2 test fixture\n' > "${A}/versions.yml"
AGG_ARGS=(--backend woltka
    --functions wtest.woltka_functions.tsv wuniq.woltka_functions.tsv
    --summaries wtest.woltka_summary.tsv wuniq.woltka_summary.tsv
    --unassigned wtest.woltka_unassigned.tsv wuniq.woltka_unassigned.tsv)
expect_ok "aggregate two samples, AGS for one" "${A}/run.log" \
    in_report "${A}" "${AGG_ARGS[@]}" --ags wtest.ags.tsv --versions-yml versions.yml --output-dir out
O="${A}/out"
check "5.3 header" \
    "[ \"\$(head -1 ${O}/read_function_abundance.tsv)\" = \"\$(printf 'sample_id\tsource\tbackend\tontology\taccession\tdescription\tabundance_tpm\tabundance_cpge\tabundance_native\tnative_unit')\" ]"
check "every row source=reads, backend=woltka, native_unit=reads, empty tpm" \
    "[ \"\$(awk -F'\t' 'NR>1 && !(\$2==\"reads\" && \$3==\"woltka\" && \$10==\"reads\" && \$7==\"\")' ${O}/read_function_abundance.tsv | wc -l)\" -eq 0 ]"
check "abundance_native = read count (wtest K00001 6.000000)" \
    "[ \"\$(value ${O}/read_function_abundance.tsv '\$1==\"wtest\" && \$5==\"K00001\"' 9)\" = '6.000000' ]"
check "abundance_cpge = RPK / GE exactly (K00003: 6.666667 / 2 = 3.333333)" \
    "[ \"\$(value ${O}/read_function_abundance.tsv '\$1==\"wtest\" && \$5==\"K00003\"' 8)\" = '3.333333' ]"
check "AGS-less sample: abundance_cpge blank, never 0" \
    "[ \"\$(awk -F'\t' '\$1==\"wuniq\" && \$8!=\"\"' ${O}/read_function_abundance.tsv | wc -l)\" -eq 0 ]"
check "all ten wide matrices written" "[ \"\$(ls ${O}/read_function_wide_*.tsv | wc -l)\" -eq 10 ]"
check "wide native: K00001 6.000000 / 0.000000 (absent under --uniq)" \
    "[ \"\$(value ${O}/read_function_wide_ko_native.tsv '\$1==\"K00001\"' 2)\" = '6.000000' ] && [ \"\$(value ${O}/read_function_wide_ko_native.tsv '\$1==\"K00001\"' 3)\" = '0.000000' ]"
check "wide cpge: AGS-less column blank" \
    "[ -z \"\$(value ${O}/read_function_wide_ko_cpge.tsv 'NR>1' 3 | tr -d '\n')\" ]"
check "sample summary: wtest ok, wuniq unavailable with 6 ambiguous reads" \
    "[ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"wtest\"' 3)\" = '0.1.7' ] && [ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"wtest\"' 9)\" = 'ok' ] && [ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"wuniq\"' 6)\" = '6.000000' ] && [ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"wuniq\"' 9)\" = 'unavailable' ]"
check "annotated fraction: wtest any 0.9000" \
    "[ \"\$(value ${O}/read_annotated_fraction.tsv '\$1==\"wtest\" && \$3==\"any\"' 6)\" = '0.9000' ]"

NA="${TMP_DIR}/agg_noags"; mkdir -p "${NA}/out"; cp "${A}"/*.tsv "${A}/versions.yml" "${NA}/"
expect_ok "no --ags at all: still succeeds" "${NA}/run.log" \
    in_report "${NA}" "${AGG_ARGS[@]}" --versions-yml versions.yml --output-dir out
check "no --ags: every cpge_status unavailable" \
    "[ \"\$(grep -c 'unavailable\$' ${NA}/out/read_sample_summary.tsv)\" -eq 2 ]"

G="${TMP_DIR}/agg_neg"; mkdir -p "${G}/out"; cp "${A}"/*.tsv "${A}/versions.yml" "${G}/"
sed 's/0.1.7/0.1.8/' "${G}/versions.yml" > "${G}/v_unknown.yml"
grep -v 'woltka:' "${G}/versions.yml" > "${G}/v_missing.yml"
sed '1s/rpk/RPK/' "${G}/wtest.woltka_functions.tsv" > "${G}/bad.woltka_functions.tsv"
mkdir -p "${G}/badont"; awk -F'\t' 'BEGIN{OFS="\t"} NR==2{$2="go"} {print}' "${G}/wtest.woltka_functions.tsv" > "${G}/badont/wtest.woltka_functions.tsv"
printf 'sample_id\taverage_genome_size_bp\tgenome_equivalents\ttotal_bases\nzzz\t3000000\t2.0\t6000000\n' > "${G}/zzz.ags.tsv"
expect_err "guard: unknown woltka version fails" "${G}/v1.log" "not among the versions" \
    in_report "${G}" "${AGG_ARGS[@]}" --versions-yml v_unknown.yml --output-dir out
expect_err "guard: missing woltka version fails" "${G}/v2.log" "no 'woltka:' version" \
    in_report "${G}" "${AGG_ARGS[@]}" --versions-yml v_missing.yml --output-dir out
expect_err "guard: renamed functions header fails" "${G}/h.log" "unexpected header" \
    in_report "${G}" --backend woltka --functions bad.woltka_functions.tsv --summaries wtest.woltka_summary.tsv \
        --unassigned wtest.woltka_unassigned.tsv --versions-yml versions.yml --output-dir out
expect_err "guard: unknown ontology fails" "${G}/o.log" "unknown ontology" \
    in_report "${G}" --backend woltka --functions badont/wtest.woltka_functions.tsv --summaries wtest.woltka_summary.tsv \
        --unassigned wtest.woltka_unassigned.tsv --versions-yml versions.yml --output-dir out
expect_err "guard: sample set mismatch fails" "${G}/s.log" "every sample needs all three" \
    in_report "${G}" --backend woltka --functions wtest.woltka_functions.tsv wuniq.woltka_functions.tsv \
        --summaries wtest.woltka_summary.tsv --unassigned wtest.woltka_unassigned.tsv wuniq.woltka_unassigned.tsv \
        --versions-yml versions.yml --output-dir out
expect_err "guard: ags for an unknown sample fails" "${G}/a.log" "unknown sample" \
    in_report "${G}" "${AGG_ARGS[@]}" --ags zzz.ags.tsv --versions-yml versions.yml --output-dir out
expect_err "guard: duplicate sample input fails" "${G}/d.log" "duplicate functions input" \
    in_report "${G}" --backend woltka --functions wtest.woltka_functions.tsv wtest.woltka_functions.tsv \
        --summaries wtest.woltka_summary.tsv --unassigned wtest.woltka_unassigned.tsv --versions-yml versions.yml --output-dir out

#
# aggregate_read_functions.py, superfocus backend (T8b, Q17)
#
echo "--- aggregate_read_functions.py: superfocus backend ---"
SFIX="${REPO_DIR}/tests/data/superfocus"
SF="${TMP_DIR}/agg_sf"; mkdir -p "${SF}/out"
cp "${SFIX}"/sftest.superfocus_{functions,summary}.tsv "${SFIX}/sftest.ags.tsv" "${SFIX}/superfocus_versions.yml" "${SF}/"
for f in "${SF}"/sftest.superfocus_*.tsv; do
    sed 's/^sftest\t/sfcopy\t/' "${f}" > "${SF}/sfcopy.${f##*/sftest.}"
done
SF_ARGS=(--backend superfocus
    --functions sftest.superfocus_functions.tsv sfcopy.superfocus_functions.tsv
    --summaries sftest.superfocus_summary.tsv sfcopy.superfocus_summary.tsv)
expect_ok "superfocus: aggregate two samples, AGS for one, no --unassigned" "${SF}/run.log" \
    in_report "${SF}" "${SF_ARGS[@]}" --ags sftest.ags.tsv --versions-yml superfocus_versions.yml --output-dir out
O="${SF}/out"
check "superfocus: every row source=reads, backend=superfocus, native_unit=reads, empty tpm AND cpge" \
    "[ \"\$(awk -F'\t' 'NR>1 && !(\$2==\"reads\" && \$3==\"superfocus\" && \$10==\"reads\" && \$7==\"\" && \$8==\"\")' ${O}/read_function_abundance.tsv | wc -l)\" -eq 0 ]"
check "superfocus: only seed_level1/2/3 ontologies (absent from KO/EC tables, 10.2)" \
    "[ \"\$(awk -F'\t' 'NR>1 {print \$4}' ${O}/read_function_abundance.tsv | sort -u | tr '\n' ' ')\" = 'seed_level1 seed_level2 seed_level3 ' ]"
check "superfocus: abundance_native = composed count (TCA cycle in plants 7.000000)" \
    "[ \"\$(value ${O}/read_function_abundance.tsv '\$1==\"sftest\" && \$6==\"TCA cycle in plants\"' 9)\" = '7.000000' ]"
check "superfocus: exactly the three native wide matrices, no _cpge, no ko/ec" \
    "[ \"\$(cd ${O} && ls read_function_wide_*.tsv | tr '\n' ' ')\" = 'read_function_wide_seed_level1_native.tsv read_function_wide_seed_level2_native.tsv read_function_wide_seed_level3_native.tsv ' ]"
check "superfocus: wide level-2 matrix keeps the '-' rows apart per level-1 parent" \
    "[ \"\$(value ${O}/read_function_wide_seed_level2_native.tsv '\$1==\"Carbohydrates | -\"' 2)\" = '2.000000' ]"
check "superfocus: sample summary not_applicable, AGS still reported, ambiguity column blank" \
    "[ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"sftest\"' 9)\" = 'not_applicable' ] && [ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"sfcopy\"' 9)\" = 'not_applicable' ] && [ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"sftest\"' 8)\" = '2.0' ] && [ -z \"\$(value ${O}/read_sample_summary.tsv '\$1==\"sftest\"' 6)\" ] && [ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"sftest\"' 3)\" = '1.8' ]"
check "superfocus: annotated fraction any 0.8636 (19 of 22 reads)" \
    "[ \"\$(value ${O}/read_annotated_fraction.tsv '\$1==\"sftest\" && \$3==\"any\"' 6)\" = '0.8636' ]"
check "superfocus: no CPGE warning printed (not applicable, not a fallback)" \
    "! grep -q 'abundance_cpge is left blank' ${SF}/run.log"

SG="${TMP_DIR}/agg_sf_neg"; mkdir -p "${SG}/out"; cp "${SF}"/*.tsv "${SF}/superfocus_versions.yml" "${SG}/"
sed 's/superfocus: 1.8/superfocus: 1.6/' "${SG}/superfocus_versions.yml" > "${SG}/v_old.yml"
grep -v 'superfocus: ' "${SG}/superfocus_versions.yml" > "${SG}/v_missing.yml"
mkdir -p "${SG}/rpk" "${SG}/ko"
awk -F'\t' 'BEGIN{OFS="\t"} NR==2{$6="1.5"} {print}' "${SG}/sftest.superfocus_functions.tsv" > "${SG}/rpk/sftest.superfocus_functions.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} NR==2{$2="ko"} {print}' "${SG}/sftest.superfocus_functions.tsv" > "${SG}/ko/sftest.superfocus_functions.tsv"
printf 'sample_id\treads_unassigned\nsftest\t0\n' > "${SG}/sftest.superfocus_unassigned.tsv"
cp "${A}/versions.yml" "${SG}/woltka_versions.yml"
expect_err "superfocus guard: unverified superfocus version fails" "${SG}/v1.log" "KNOWN_SUPERFOCUS_VERSIONS" \
    in_report "${SG}" "${SF_ARGS[@]}" --versions-yml v_old.yml --output-dir out
expect_err "superfocus guard: missing superfocus version fails" "${SG}/v2.log" "no 'superfocus:' version" \
    in_report "${SG}" "${SF_ARGS[@]}" --versions-yml v_missing.yml --output-dir out
expect_err "superfocus guard: woltka versions.yml for a superfocus run fails" "${SG}/v3.log" "no 'superfocus:' version" \
    in_report "${SG}" "${SF_ARGS[@]}" --versions-yml woltka_versions.yml --output-dir out
expect_err "superfocus guard: filled rpk fails (CPGE not computable)" "${SG}/r.log" "rpk must be blank" \
    in_report "${SG}" --backend superfocus --functions rpk/sftest.superfocus_functions.tsv \
        --summaries sftest.superfocus_summary.tsv --versions-yml superfocus_versions.yml --output-dir out
expect_err "superfocus guard: a KO row in a superfocus table fails (never mapped to KO)" "${SG}/k.log" "unknown ontology" \
    in_report "${SG}" --backend superfocus --functions ko/sftest.superfocus_functions.tsv \
        --summaries sftest.superfocus_summary.tsv --versions-yml superfocus_versions.yml --output-dir out
expect_err "superfocus guard: --unassigned is rejected for superfocus" "${SG}/u.log" "not used by the superfocus backend" \
    in_report "${SG}" --backend superfocus --functions sftest.superfocus_functions.tsv \
        --summaries sftest.superfocus_summary.tsv --unassigned sftest.superfocus_unassigned.tsv \
        --versions-yml superfocus_versions.yml --output-dir out
expect_err "woltka guard: missing --unassigned fails" "${SG}/w.log" "are required for the woltka backend" \
    in_report "${A}" --backend woltka --functions wtest.woltka_functions.tsv --summaries wtest.woltka_summary.tsv \
        --versions-yml versions.yml --output-dir out

echo "--- aggregate_read_functions.py: humann backend ---"
HFIX="${REPO_DIR}/tests/data/humann"
HU="${TMP_DIR}/agg_hu"; mkdir -p "${HU}/out"
cp "${HFIX}"/htest.humann_{functions,summary}.tsv "${HFIX}/htest.ags.tsv" "${HFIX}/humann_versions.yml" "${HU}/"
for f in "${HU}"/htest.humann_*.tsv; do
    sed 's/^htest\t/hcopy\t/' "${f}" > "${HU}/hcopy.${f##*/htest.}"
done
HU_ARGS=(--backend humann
    --functions htest.humann_functions.tsv hcopy.humann_functions.tsv
    --summaries htest.humann_summary.tsv hcopy.humann_summary.tsv)
expect_ok "humann: aggregate two samples, AGS for one, no --unassigned" "${HU}/run.log" \
    in_report "${HU}" "${HU_ARGS[@]}" --ags htest.ags.tsv --versions-yml humann_versions.yml --output-dir out
O="${HU}/out"
check "humann: every row source=reads, backend=humann, native_unit=rpk, empty tpm" \
    "[ \"\$(awk -F'\t' 'NR>1 && !(\$2==\"reads\" && \$3==\"humann\" && \$10==\"rpk\" && \$7==\"\")' ${O}/read_function_abundance.tsv | wc -l)\" -eq 0 ]"
check "humann: only ko, ec and metacyc ontologies" \
    "[ \"\$(awk -F'\t' 'NR>1 {print \$4}' ${O}/read_function_abundance.tsv | sort -u | tr '\n' ' ')\" = 'ec ko metacyc ' ]"
K2=$(value "${HFIX}/htest.humann_functions.tsv" '$2=="ko" && $3=="K00002"' 6)
check "humann: abundance_native = RPK (K00002)" \
    "[ \"\$(value ${O}/read_function_abundance.tsv '\$1==\"htest\" && \$5==\"K00002\"' 9)\" = \"\$(printf '%.6f' ${K2})\" ]"
check "humann: abundance_cpge = RPK / GE exactly (GE 2.0)" \
    "[ \"\$(value ${O}/read_function_abundance.tsv '\$1==\"htest\" && \$5==\"K00002\"' 8)\" = \"\$(python3 -c \"print('%.6f' % (${K2} / 2.0))\")\" ]"
check "humann: AGS-less sample has blank CPGE, never 0" \
    "[ \"\$(awk -F'\t' '\$1==\"hcopy\" && \$8!=\"\"' ${O}/read_function_abundance.tsv | wc -l)\" -eq 0 ]"
check "humann: six wide matrices (ko/ec/metacyc x native/cpge)" \
    "[ \"\$(cd ${O} && ls read_function_wide_*.tsv | tr '\n' ' ')\" = 'read_function_wide_ec_cpge.tsv read_function_wide_ec_native.tsv read_function_wide_ko_cpge.tsv read_function_wide_ko_native.tsv read_function_wide_metacyc_cpge.tsv read_function_wide_metacyc_native.tsv ' ]"
check "humann: annotated fraction 'any' from reads, ontology rows = composer RPK share with blank read columns" \
    "[ \"\$(value ${O}/read_annotated_fraction.tsv '\$1==\"htest\" && \$3==\"any\"' 6)\" = \"\$(printf '%.4f' $(value "${HFIX}/htest.humann_summary.tsv" '$2=="any"' 5))\" ] && [ \"\$(value ${O}/read_annotated_fraction.tsv '\$1==\"htest\" && \$3==\"ko\"' 6)\" = \"\$(printf '%.4f' $(value "${HFIX}/htest.humann_summary.tsv" '$2=="ko"' 5))\" ] && [ -z \"\$(value ${O}/read_annotated_fraction.tsv '\$3!=\"any\" && NR>1' 4 | tr -d '\n')\" ]"
check "humann: sample summary ok / unavailable, tool 4.0.0.alpha.2, db provenance" \
    "[ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"htest\"' 9)\" = 'ok' ] && [ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"hcopy\"' 9)\" = 'unavailable' ] && [ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"htest\"' 3)\" = '4.0.0.alpha.2' ] && [ \"\$(value ${O}/read_sample_summary.tsv '\$1==\"htest\"' 4)\" = 'fixture' ]"

HG="${TMP_DIR}/agg_hu_neg"; mkdir -p "${HG}/out" "${HG}/reads" "${HG}/frac"; cp "${HU}"/*.tsv "${HU}/humann_versions.yml" "${HG}/"
sed 's/humann: 4.0.0.alpha.2/humann: 4.0.0.alpha.3/' "${HG}/humann_versions.yml" > "${HG}/v_new.yml"
awk -F'\t' 'BEGIN{OFS="\t"} $2=="ko"{$3="12"} {print}' "${HG}/htest.humann_summary.tsv" > "${HG}/reads/htest.humann_summary.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} $2=="ko"{$5="1.5"} {print}' "${HG}/htest.humann_summary.tsv" > "${HG}/frac/htest.humann_summary.tsv"
expect_err "humann guard: unverified HUMAnN version fails" "${HG}/v.log" "KNOWN_HUMANN_VERSIONS" \
    in_report "${HG}" "${HU_ARGS[@]}" --versions-yml v_new.yml --output-dir out
expect_err "humann guard: read counts on an ontology summary row fail" "${HG}/r.log" "must have blank read columns" \
    in_report "${HG}" --backend humann --functions htest.humann_functions.tsv \
        --summaries reads/htest.humann_summary.tsv --versions-yml humann_versions.yml --output-dir out
expect_err "humann guard: RPK share outside [0,1] fails" "${HG}/f.log" "outside \[0, 1\]" \
    in_report "${HG}" --backend humann --functions htest.humann_functions.tsv \
        --summaries frac/htest.humann_summary.tsv --versions-yml humann_versions.yml --output-dir out

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
