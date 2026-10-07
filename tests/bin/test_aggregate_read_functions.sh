#!/bin/bash
#
# Direct tests for the read-level functional branch scripts (design doc
# Sections 4.6, 4.6.2, 5.3, 6.2, 10.2; task T8a - Woltka backend), each run
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
#     - --uniq: ambiguous mates move to the unassigned table
#     - header-only (empty-alignment) profile -> header-only functions,
#       zero summary rows with empty fractions
#     - guards fail loudly: a demultiplexed (3-column) profile, a profile ORF
#       missing from length.map, a missing database file
#   bin/aggregate_read_functions.py  (pandas image, AGGREGATE_READ_FUNCTIONS)
#     - Section 5.3 schema with source=reads, backend=woltka, empty
#       abundance_tpm, abundance_native = read counts, native_unit=reads
#     - CPGE = RPK / genome_equivalents exactly; blank (never 0) for samples
#       without AGS, with cpge_status 'unavailable' in read_sample_summary
#     - wide native/cpge matrices for all five ontologies, always written
#     - version-aware guard (unknown or missing woltka version), header,
#       sample-set, ontology, duplicate and ags guards all fail loudly
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
PROFILE='python3 /pipeline_bin/woltka_function_profile.py --db /fix/db'

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
check "repeated Pfam domain counts once (PF00002.1 = 3, not 9)" \
    "[ \"\$(value ${F} '\$2==\"pfam\" && \$3==\"PF00002.1\"' 5)\" = '3' ]"
check "multi-path MetaCyc pathway counts each ORF once (PWY-1 = 7, not 10)" \
    "[ \"\$(value ${F} '\$2==\"metacyc\" && \$3==\"PWY-1\"' 5)\" = '7' ]"
check "COG ids via ko-to-cog (COG0001 = 14)" \
    "[ \"\$(value ${F} '\$2==\"cog\" && \$3==\"COG0001\"' 5)\" = '14' ]"
check "RPK = count / kb (K00003: 8 reads / 1.2 kb = 6.6667)" \
    "value ${F} '\$2==\"ko\" && \$3==\"K00003\"' 6 | grep -q '^6\.66666666'"
check "ec and cog carry no description" \
    "[ -z \"\$(value ${F} '\$2==\"ec\" || \$2==\"cog\"' 4 | tr -d '\n')\" ]"
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
    in_woltka "${N}" "mkdir -p partial && cp -r /fix/db/proteins partial/ && python3 /pipeline_bin/woltka_function_profile.py --db partial --profile orf.tsv --sample-id w --prefix w"

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

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
