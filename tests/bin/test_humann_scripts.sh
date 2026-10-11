#!/bin/bash
#
# Direct tests for the HUMAnN read-level functional backend (design doc
# Sections 4.6.1, 5.3, 10.2; task T8c, Q2), run inside the exact image the
# HUMANN module pins:
#
#   real `humann` (4.0.0a2) on the committed fixture (tests/data/humann, made
#   by make_humann_fixtures.sh) with the module's database flags plus the
#   fixture arguments (--bypass-prescreen, subject-coverage thresholds 0):
#     - the raw tables match the committed ones byte for byte
#   bin/humann_function_profile.py
#     - KO / EC regrouping with the synthetic maps: a family carrying two KOs
#       counts toward both, once each; a KO fed by three families sums them;
#       UniClust90 families reach EC, never KO; READS_UNMAPPED never reaches a
#       term (the HUMAnN 4.0.0a2 regroup bug the composer avoids)
#     - metacyc rows = unstratified pathways minus UNMAPPED / UNINTEGRATED
#     - summary: 'any' reads (total, total - READS_UNMAPPED); per-ontology RPK
#       shares with blank read columns
#     - no-reads mode (no tables): header-only functions, zero summary
#     - version-aware guards fail loudly (design doc 10.2: "deliberately
#       feeding it a malformed file"): CPM header (run without RPKs), unknown
#       HUMAnN version, wrong sample column, 3-column row, non-numeric and
#       negative values, unknown feature forms, missing map file,
#       READS_UNMAPPED above the input reads, one table without the other
#
# Requirements: docker (present on GitHub ubuntu-latest runners).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
FIX="${REPO_DIR}/tests/data/humann"

# Exact image pinned in the module, so the test cannot drift from it
HUMANN_IMG=$(grep -o "'community.wave.seqera.io/library/python_metaphlan_diamond_bowtie2[^']*'" \
    "${REPO_DIR}/modules/local/humann/main.nf" | grep -v oras | tr -d "'" | head -1)
if [ -z "${HUMANN_IMG}" ]; then
    echo "ERROR: could not extract the pinned HUMAnN image from the module"
    exit 1
fi

echo "=== HUMAnN backend script test suite ==="
echo "HUMAnN image: ${HUMANN_IMG}"
echo ""

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0

in_humann() {
    # in_humann <workdir> <shell command> - fixture mounted read-only at /fix
    local workdir="$1"; shift
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" -v "${FIX}:/fix:ro" \
        -v "${workdir}:/w" -w /w "${HUMANN_IMG}" bash -c "$*"
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
# sum of the unstratified gene-family RPK of the families listed in a map line
family_sum() {  # family_sum <genefamilies.tsv> <family...>
    local gf="$1"; shift
    awk -F'\t' -v ids="$*" 'BEGIN { n = split(ids, a, " "); for (i = 1; i <= n; i++) want[a[i]] = 1 }
        ($1 in want) { s += $2 } END { printf "%.6f", s }' "${gf}"
}

COMPOSE='python3 /pipeline_bin/humann_function_profile.py --utility-db /fix/db/utility_mapping --basename htest --sample-id htest --prefix htest'
READS=$(( $(zcat "${FIX}/reads/htest_R1.fastq.gz" "${FIX}/reads/htest_R2.fastq.gz" | wc -l) / 4 ))

#
# Real HUMAnN on the fixture
#
echo "--- humann 4.0.0a2 on the fixture (humann image) ---"
H="${TMP_DIR}/humann"; mkdir -p "${H}"
expect_ok "humann runs on the fixture with the module's database flags" "${H}/run.log" \
    in_humann "${H}" "zcat /fix/reads/htest_R1.fastq.gz /fix/reads/htest_R2.fastq.gz > reads.fastq && \
        humann --input reads.fastq --input-format fastq --output . --output-basename htest --threads 1 \
            --nucleotide-database /fix/db/chocophlan --protein-database /fix/db/uniref \
            --utility-database /fix/db/utility_mapping \
            --pathways-database /fix/db/utility_mapping/metacyc_reactions_level4ec_only.uniref.bz2,/fix/db/utility_mapping/metacyc_pathways_structured_filtered_v24_subreactions \
            --count-normalization RPKs --remove-temp-output \
            --bypass-prescreen --nucleotide-subject-coverage-threshold 0 --translated-subject-coverage-threshold 0"
for t in 2_genefamilies 3_reactions 4_pathabundance; do
    check "raw htest_${t}.tsv matches the committed fixture byte for byte" \
        "cmp -s ${H}/htest_${t}.tsv ${FIX}/htest_${t}.tsv"
done
check "gene-family header is the verified RPK layout" \
    "[ \"\$(head -1 ${H}/htest_2_genefamilies.tsv)\" = \"\$(printf '# Gene Family HUMAnN v4.0.0.alpha.2 RPKs\thtest')\" ]"
check "no pathway coverage table (writer disabled in 4.0.0a2)" "[ ! -e ${H}/htest_5_pathcoverage.tsv ]"

# The module's --metaphlan-options must be accepted by the pinned MetaPhlAn
# (4.1.2 names its database option --bowtie2db; --db_dir only exists from 4.2
# and was rejected on the real-database acceptance run). No MetaPhlAn database
# exists in CI, so the run is expected to fail on the empty placeholder files -
# but never on argument parsing. The option string is taken from the module.
MPA_OPTS=$(grep -o -- '--metaphlan-options "[^"]*"' "${REPO_DIR}/modules/local/humann/main.nf" \
    | sed -e 's/^--metaphlan-options "//' -e 's/"$//' \
          -e 's/\\\${METAPHLAN_DB_DIR}/\/fix\/db\/metaphlan/' \
          -e 's/\${mpa_index}/mpa_vOct22_CHOCOPhlAnSGB_202403/' -e 's/\${task.cpus}/1/')
M="${TMP_DIR}/mpa"; mkdir -p "${M}"
in_humann "${M}" "zcat /fix/reads/htest_R1.fastq.gz > r.fastq; metaphlan r.fastq --input_type fastq ${MPA_OPTS} --offline \
    -o profile.tsv --bowtie2out bt2.txt" > "${M}/run.log" 2>&1 || true
check "module MetaPhlAn options extracted (${MPA_OPTS})" "[ -n \"${MPA_OPTS}\" ] && echo '${MPA_OPTS}' | grep -q -- '--bowtie2db /fix/db/metaphlan --index mpa_vOct22_CHOCOPhlAnSGB_202403'"
check "MetaPhlAn 4.1.2 accepts every option the module passes" \
    "! grep -qE 'unrecognized arguments|invalid choice' ${M}/run.log"

#
# humann_function_profile.py
#
echo "--- humann_function_profile.py (humann image) ---"
D="${TMP_DIR}/compose"; mkdir -p "${D}"
expect_ok "compose the committed raw tables" "${D}/run.log" \
    in_humann "${D}" "${COMPOSE} --genefamilies /fix/htest_2_genefamilies.tsv --pathabundance /fix/htest_4_pathabundance.tsv --input-reads ${READS}"
F="${D}/htest.humann_functions.tsv"
SUM="${D}/htest.humann_summary.tsv"
GF="${FIX}/htest_2_genefamilies.tsv"
KOMAP=$(zcat "${FIX}/db/utility_mapping/map_ko_uniref90.txt.gz")
ECMAP=$(zcat "${FIX}/db/utility_mapping/map_level4ec_uniclust90.txt.gz")
members() { awk -F'\t' -v t="$2" '$1 == t { for (i = 2; i <= NF; i++) print $i }' <<< "$1" | tr '\n' ' '; }
check "functions and summary match the committed fixture outputs" \
    "cmp -s ${F} ${FIX}/htest.humann_functions.tsv && cmp -s ${SUM} ${FIX}/htest.humann_summary.tsv"
check "header matches the composer layout" \
    "[ \"\$(head -1 ${F})\" = \"\$(printf 'sample_id\tontology\taccession\tdescription\tcount\trpk')\" ]"
check "only ko, ec and metacyc rows" \
    "[ \"\$(tail -n +2 ${F} | cut -f2 | sort -u | tr '\n' ' ')\" = 'ec ko metacyc ' ]"
check "K00002 = sum of its three families' RPK" \
    "[ \"\$(printf '%.6f' \$(value ${F} '\$2==\"ko\" && \$3==\"K00002\"' 5))\" = \"\$(family_sum ${GF} \$(members \"\${KOMAP}\" K00002))\" ]"
check "a family with two KOs counts toward both, once (K00001 = its one present family)" \
    "[ \"\$(printf '%.6f' \$(value ${F} '\$2==\"ko\" && \$3==\"K00001\"' 5))\" = \"\$(family_sum ${GF} \$(members \"\${KOMAP}\" K00001))\" ]"
check "EC 2.7.7.7 sums its UniClust90 and UniRef90 families" \
    "[ \"\$(printf '%.6f' \$(value ${F} '\$2==\"ec\" && \$3==\"2.7.7.7\"' 5))\" = \"\$(family_sum ${GF} \$(members \"\${ECMAP}\" 2.7.7.7))\" ]"
check "count = rpk on every row (HUMAnN RPK is the native unit)" \
    "[ -z \"\$(awk -F'\t' 'NR>1 && \$5 != \$6' ${F})\" ]"
check "names attached from the name maps" \
    "[ \"\$(value ${F} '\$2==\"ko\" && \$3==\"K00003\"' 4)\" = 'homoserine dehydrogenase [EC:1.1.1.3]' ]"
check "metacyc rows = unstratified pathways minus UNMAPPED/UNINTEGRATED" \
    "[ \"\$(value ${F} '\$2==\"metacyc\"' 3 | wc -l)\" -eq \"\$(grep -v '|' ${FIX}/htest_4_pathabundance.tsv | tail -n +2 | grep -cvE '^(UNMAPPED|UNINTEGRATED)')\" ]"
UNMAPPED=$(awk -F'\t' '$1=="READS_UNMAPPED" { printf "%d", $2 }' "${GF}")
check "summary any: reads given / reads minus READS_UNMAPPED" \
    "[ \"\$(value ${SUM} '\$2==\"any\"' 3)\" = '${READS}' ] && [ \"\$(value ${SUM} '\$2==\"any\"' 4)\" = \"\$((READS - UNMAPPED))\" ]"
check "summary ontology rows: blank read columns, RPK share in [0,1]" \
    "[ -z \"\$(awk -F'\t' 'NR>2 && (\$3 != \"\" || \$4 != \"\" || \$5 < 0 || \$5 > 1)' ${SUM})\" ]"

E="${TMP_DIR}/empty"; mkdir -p "${E}"
expect_ok "no-reads mode (no tables)" "${E}/run.log" in_humann "${E}" "${COMPOSE} --input-reads 0"
check "no reads: header-only functions" "[ \"\$(wc -l < ${E}/htest.humann_functions.tsv)\" -eq 1 ]"
check "no reads: zero 'any' row, blank fractions" \
    "[ \"\$(value ${E}/htest.humann_summary.tsv '\$2==\"any\"' 3)\" = '0' ] && [ -z \"\$(value ${E}/htest.humann_summary.tsv 'NR>1' 5 | tr -d '\n')\" ]"

#
# Version-aware guards (design doc 10.2): deliberately malformed inputs
#
echo "--- humann_function_profile.py guards ---"
N="${TMP_DIR}/neg"; mkdir -p "${N}"
cp "${FIX}/htest_2_genefamilies.tsv" "${FIX}/htest_4_pathabundance.tsv" "${N}/"
sed '1s/RPKs/Adjusted CPMs/' "${N}/htest_2_genefamilies.tsv" > "${N}/cpm.tsv"
sed '1s/v4.0.0.alpha.2/v4.0.0.alpha.3/' "${N}/htest_2_genefamilies.tsv" > "${N}/version.tsv"
sed '1s/\thtest$/\tother/' "${N}/htest_2_genefamilies.tsv" > "${N}/column.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} NR==3 { $0 = $0 "\textra" } { print }' "${N}/htest_2_genefamilies.tsv" > "${N}/cols.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} NR==3 { $2 = "abc" } { print }' "${N}/htest_2_genefamilies.tsv" > "${N}/nonnum.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} NR==3 { $2 = "-1" } { print }' "${N}/htest_2_genefamilies.tsv" > "${N}/neg.tsv"
# inserted mid-file: HUMAnN writes no trailing newline, so appending would merge rows
awk 'NR==3 { print "UniProt_P12345\t1.0" } { print }' "${N}/htest_2_genefamilies.tsv" > "${N}/feature.tsv"
awk 'NR==3 { print "MYSTERY\t1.0" } { print }' "${N}/htest_4_pathabundance.tsv" > "${N}/pwyfeature.tsv"
sed '1s/htest_Abundance/htest_Coverage/' "${N}/htest_4_pathabundance.tsv" > "${N}/pwycolumn.tsv"
mkdir -p "${N}/nomap" && cp "${FIX}"/db/utility_mapping/map_level4ec_* "${FIX}"/db/utility_mapping/map_ko_name.txt.gz "${N}/nomap/"
G='--pathabundance /w/htest_4_pathabundance.tsv --input-reads '"${READS}"
expect_err "guard: CPM table (run without --count-normalization RPKs) fails" "${N}/cpm.log" "unexpected gene-family header" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/cpm.tsv ${G}"
expect_err "guard: unverified HUMAnN version fails" "${N}/version.log" "not among the versions this parser was verified against" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/version.tsv ${G}"
expect_err "guard: wrong sample column fails" "${N}/column.log" "sample column" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/column.tsv ${G}"
expect_err "guard: 3-column row fails" "${N}/cols.log" "expected 2 columns" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/cols.tsv ${G}"
expect_err "guard: non-numeric value fails" "${N}/nonnum.log" "non-numeric value" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/nonnum.tsv ${G}"
expect_err "guard: negative value fails" "${N}/neg.log" "negative value" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/neg.tsv ${G}"
expect_err "guard: unknown gene-family feature form fails" "${N}/feature.log" "unexpected gene-family feature" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/feature.tsv ${G}"
expect_err "guard: unknown pathway feature form fails" "${N}/pwyfeature.log" "unexpected pathway feature" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/htest_2_genefamilies.tsv --pathabundance /w/pwyfeature.tsv --input-reads ${READS}"
expect_err "guard: wrong pathway sample column fails" "${N}/pwycolumn.log" "sample column" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/htest_2_genefamilies.tsv --pathabundance /w/pwycolumn.tsv --input-reads ${READS}"
expect_err "guard: READS_UNMAPPED above the input reads fails" "${N}/reads.log" "exceeds the" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/htest_2_genefamilies.tsv ${G%${READS}}1"
expect_err "guard: missing KO map fails" "${N}/map.log" "utility mapping file missing" \
    in_humann "${N}" "python3 /pipeline_bin/humann_function_profile.py --utility-db /w/nomap --basename htest --sample-id htest --prefix htest --genefamilies /w/htest_2_genefamilies.tsv ${G}"
expect_err "guard: one table without the other fails" "${N}/half.log" "must be given together" \
    in_humann "${N}" "${COMPOSE} --genefamilies /w/htest_2_genefamilies.tsv --input-reads ${READS}"

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
