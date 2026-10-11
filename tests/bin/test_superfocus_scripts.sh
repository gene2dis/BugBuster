#!/bin/bash
#
# Direct tests for the SUPER-FOCUS read-level functional backend (design doc
# Sections 4.6.3, 5.3, 10.2; task T8b, Q17), run inside the exact image the
# SUPERFOCUS module pins:
#
#   real `superfocus` (1.8) on the committed fixture (tests/data/superfocus,
#   made by make_superfocus_fixtures.sh), invoked as the module does, with
#   BOTH aligners (diamond, mmseqs2):
#     - concatenated R1+R2+S query: mates count as separate reads
#       (22 input reads, 19 with an accepted hit)
#     - the fumarate hydratase reads tie between two subsystems and are
#       divided 1/k (3 + 3 of 6 reads)
#     - SUPER-FOCUS's own level-2 file merges the '-' placeholder across
#       level-1 parents (the reason the pipeline recomputes levels)
#     - an all-unrelated query gives a header-only table (exit 0)
#   bin/superfocus_function_profile.py
#     - per-level sums with path-qualified accessions, '-' level 2 kept
#       apart per level-1 parent, rpk blank, summary fractions
#     - no-reads mode (no --table): header-only functions, zero summary
#     - layout guards fail loudly: wrong database / aligner line, separator
#       row, header (query column name), field count, non-numeric /
#       negative count, empty level, a total above the input reads (-n 0)
#
# Requirements: docker (present on GitHub ubuntu-latest runners).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
FIX="${REPO_DIR}/tests/data/superfocus"

# Exact image pinned in the module, so the test cannot drift from it
SF_IMG=$(grep -o "'community.wave.seqera.io/library/super-focus[^']*'" "${REPO_DIR}/modules/local/superfocus/main.nf" \
    | tr -d "'" | head -1)
if [ -z "${SF_IMG}" ]; then
    echo "ERROR: could not extract the pinned SUPER-FOCUS image from the module"
    exit 1
fi

echo "=== SUPER-FOCUS backend script test suite ==="
echo "SUPER-FOCUS image: ${SF_IMG}"
echo ""

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0

in_sf() {
    # in_sf <workdir> <shell command> - fixture mounted read-only at /fix
    local workdir="$1"; shift
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" -v "${FIX}:/fix:ro" \
        -v "${workdir}:/w" -w /w "${SF_IMG}" bash -c "$*"
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

COMPOSE='python3 /pipeline_bin/superfocus_function_profile.py'

#
# real superfocus, both aligners (module invocation)
#
echo "--- superfocus 1.8 on the fixture (both aligners) ---"
for aligner in diamond mmseqs2; do
    R="${TMP_DIR}/run_${aligner}"; mkdir -p "${R}"
    expect_ok "${aligner}: superfocus + composer succeed" "${R}/run.log" in_sf "${R}" "
        set -euo pipefail
        gzip -cdf /fix/reads/sftest_R1.fastq.gz /fix/reads/sftest_R2.fastq.gz /fix/reads/sftest_S.fastq.gz > sftest.fastq
        n=\$(( \$(wc -l < sftest.fastq) / 4 ))
        superfocus -q sftest.fastq -dir out -a ${aligner} -db DB_90 -b /fix/sf_db -t 1 -tmp ./tmp -n 1 -d -l sf.log \
            -mi 60 -ml 15 -e 0.00001 -f 1
        ${COMPOSE} --table out/output_all_levels_and_function.xls --query sftest.fastq --input-reads \${n} \
            --aligner ${aligner} --database 90 --sample-id sftest --prefix sftest"
    check "${aligner}: composed functions identical to the committed fixture" \
        "cmp -s ${R}/sftest.superfocus_functions.tsv ${FIX}/sftest.superfocus_functions.tsv"
    check "${aligner}: summary identical to the committed fixture (22 reads, 19 hit)" \
        "cmp -s ${R}/sftest.superfocus_summary.tsv ${FIX}/sftest.superfocus_summary.tsv"
    check "${aligner}: tie divided 1/k (TCA Cycle fumarate hydratase 3.0 of 6 reads)" \
        "[ \"\$(awk -F'\t' '\$3==\"TCA Cycle\" {print \$5}' ${R}/out/output_all_levels_and_function.xls)\" = '3.0' ]"
done
check "upstream level-2 file merges '-' across level-1 parents (5 = 3 + 2)" \
    "[ \"\$(awk -F'\t' '\$1==\"-\" {print \$2}' ${TMP_DIR}/run_diamond/out/output_subsystem_level_2.xls)\" = '5' ]"

N="${TMP_DIR}/nohit"; mkdir -p "${N}"
expect_ok "all-unrelated query: superfocus exits 0" "${N}/run.log" in_sf "${N}" "
    set -euo pipefail
    gzip -cdf /fix/reads/sftest_S.fastq.gz | tail -4 > nohit.fastq
    superfocus -q nohit.fastq -dir out -a diamond -db DB_90 -b /fix/sf_db -t 1 -tmp ./tmp -n 1 -d -l sf.log
    ${COMPOSE} --table out/output_all_levels_and_function.xls --query nohit.fastq --input-reads 1 \
        --aligner diamond --database 90 --sample-id nohit --prefix nohit"
check "all-unrelated query: header-only functions, 0 of 1 annotated" \
    "[ \"\$(wc -l < ${N}/nohit.superfocus_functions.tsv)\" -eq 1 ] && [ \"\$(value ${N}/nohit.superfocus_summary.tsv '\$2==\"any\"' 4)\" = '0' ] && [ \"\$(value ${N}/nohit.superfocus_summary.tsv '\$2==\"any\"' 5)\" = '0.000000' ]"

#
# composer on the committed raw table
#
echo "--- superfocus_function_profile.py ---"
C="${TMP_DIR}/compose"; mkdir -p "${C}"
cp "${FIX}/sftest.superfocus_all_levels_and_function.xls" "${C}/raw.xls"
F="${FIX}/sftest.superfocus_functions.tsv"
check "functions header" \
    "[ \"\$(head -1 ${F})\" = \"\$(printf 'sample_id\tontology\taccession\tdescription\tcount\trpk')\" ]"
check "seed_level1: Carbohydrates 16, Amino Acids and Derivatives 3" \
    "[ \"\$(value ${F} '\$2==\"seed_level1\" && \$3==\"Carbohydrates\"' 5)\" = '16' ] && [ \"\$(value ${F} '\$2==\"seed_level1\" && \$3==\"Amino Acids and Derivatives\"' 5)\" = '3' ]"
check "seed_level2: '-' kept apart per level-1 parent (Carbohydrates | - 2; Amino Acids and Derivatives | - 3)" \
    "[ \"\$(value ${F} '\$3==\"Carbohydrates | -\"' 5)\" = '2' ] && [ \"\$(value ${F} '\$3==\"Amino Acids and Derivatives | -\"' 5)\" = '3' ]"
check "seed_level3: TCA cycle in plants 7 (3 tie + 4 SDH1), description = level name" \
    "[ \"\$(value ${F} '\$3==\"Carbohydrates | Central carbohydrate metabolism | TCA cycle in plants\"' 5)\" = '7' ] && [ \"\$(value ${F} '\$3==\"Carbohydrates | Central carbohydrate metabolism | TCA cycle in plants\"' 4)\" = 'TCA cycle in plants' ]"
check "every level sums to the 19 hit reads" \
    "[ \"\$(awk -F'\t' 'NR>1 {s[\$2]+=\$5} END {print s[\"seed_level1\"], s[\"seed_level2\"], s[\"seed_level3\"]}' ${F})\" = '19 19 19' ]"
check "rpk always blank (no gene length, CPGE not applicable)" \
    "[ \"\$(awk -F'\t' 'NR>1 && \$6!=\"\"' ${F} | wc -l)\" -eq 0 ]"
check "summary: any + three levels, 19/22 = 0.863636" \
    "[ \"\$(cut -f2 ${FIX}/sftest.superfocus_summary.tsv | tail -n +2 | tr '\n' ' ')\" = 'any seed_level1 seed_level2 seed_level3 ' ] && [ \"\$(value ${FIX}/sftest.superfocus_summary.tsv '\$2==\"any\"' 5)\" = '0.863636' ]"

expect_ok "no-reads mode (no --table)" "${C}/empty.log" \
    in_sf "${C}" "${COMPOSE} --query e.fastq --input-reads 0 --aligner diamond --database 90 --sample-id e --prefix e"
check "no-reads mode: header-only functions, zero summary with empty fractions" \
    "[ \"\$(wc -l < ${C}/e.superfocus_functions.tsv)\" -eq 1 ] && [ \"\$(value ${C}/e.superfocus_summary.tsv '\$2==\"any\"' 3)\" = '0' ] && [ -z \"\$(value ${C}/e.superfocus_summary.tsv '\$2==\"any\"' 5)\" ]"

GOOD="--query sftest.fastq --input-reads 22 --aligner diamond --database 90 --sample-id g --prefix g"
mutate() {  # mutate <name> <sed/awk shell command reading raw.xls>
    in_sf "${C}" "$2 > $1.xls"
}
mutate db "sed '2s/90/95/' raw.xls"
mutate aligner "sed '3s/diamond/rapsearch/' raw.xls"
mutate sep "sed '4s/.*//' raw.xls"
mutate query "sed '5s/sftest.fastq/sftest_R1.fastq/g' raw.xls"
mutate fields "awk 'BEGIN{FS=OFS=\"\\t\"} NR==6{\$7=\"x\"} {print}' raw.xls"
mutate nonnum "awk 'BEGIN{FS=OFS=\"\\t\"} NR==6{\$5=\"three\"} {print}' raw.xls"
mutate negative "awk 'BEGIN{FS=OFS=\"\\t\"} NR==6{\$5=\"-3\"} {print}' raw.xls"
mutate emptylevel "awk 'BEGIN{FS=OFS=\"\\t\"} NR==6{\$1=\"\"} {print}' raw.xls"
expect_err "guard: database line mismatch fails" "${C}/g1.log" "expected 'Database used: 90'" \
    in_sf "${C}" "${COMPOSE} --table db.xls ${GOOD}"
expect_err "guard: aligner line mismatch fails" "${C}/g2.log" "expected 'Aligner used: diamond'" \
    in_sf "${C}" "${COMPOSE} --table aligner.xls ${GOOD}"
expect_err "guard: mmseqs2 run checked against 'mmseqs'" "${C}/g3.log" "expected 'Aligner used: mmseqs'" \
    in_sf "${C}" "${COMPOSE} --table raw.xls --query sftest.fastq --input-reads 22 --aligner mmseqs2 --database 90 --sample-id g --prefix g"
expect_err "guard: separator row changed fails" "${C}/g4.log" "blank separator row" \
    in_sf "${C}" "${COMPOSE} --table sep.xls ${GOOD}"
expect_err "guard: count column not named after the single query fails" "${C}/g5.log" "unexpected header" \
    in_sf "${C}" "${COMPOSE} --table query.xls ${GOOD}"
expect_err "guard: extra column fails" "${C}/g6.log" "fields, expected" \
    in_sf "${C}" "${COMPOSE} --table fields.xls ${GOOD}"
expect_err "guard: non-numeric count fails" "${C}/g7.log" "non-numeric count" \
    in_sf "${C}" "${COMPOSE} --table nonnum.xls ${GOOD}"
expect_err "guard: negative count fails" "${C}/g8.log" "negative count" \
    in_sf "${C}" "${COMPOSE} --table negative.xls ${GOOD}"
expect_err "guard: empty subsystem level fails" "${C}/g9.log" "empty subsystem level" \
    in_sf "${C}" "${COMPOSE} --table emptylevel.xls ${GOOD}"
expect_err "guard: total above the input reads (-n 0 / wrong query) fails" "${C}/g10.log" "exceeds the 10 input reads" \
    in_sf "${C}" "${COMPOSE} --table raw.xls --query sftest.fastq --input-reads 10 --aligner diamond --database 90 --sample-id g --prefix g"
expect_err "guard: reads but no --table fails" "${C}/g11.log" "--table is required" \
    in_sf "${C}" "${COMPOSE} --query sftest.fastq --input-reads 5 --aligner diamond --database 90 --sample-id g --prefix g"

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
