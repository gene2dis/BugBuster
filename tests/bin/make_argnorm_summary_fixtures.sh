#!/usr/bin/env bash
# Regenerate tests/data/summary/argnorm/ (roadmap Phase 27b): argNorm 1.1.0 (the image ARGNORM pins) on
#  - a hand-written ARGs-OAP 4.0.0 stage two table of sampleA (header row, 18 columns; real SARG ids, synthetic
#    reads): KPC and TEM, which the K-MARVEL / RGI fixtures of sampleA also call; mecA on 3 Strict hits (ARGs-OAP
#    only, tier 1) and 1 Loose hit of the same subtype; vanA on 2 hits (under the default --summary_args_oap_min_hits 3); rosB on 3 Loose hits; tetA(58) on
#    3 Strict hits, whose ARO card/metadata.tsv leaves out on purpose;
#  - the DeepARG tables of tests/data/summary/stages/assembly/ (sampleA, sampleB; Strict, Perfect and Loose
#    mappings, some on binned contigs);
# and card/, a hand-written CARD-metadata table for those AROs (below). Needs docker; downloads nothing but the
# argNorm image.
# usage: tests/bin/make_argnorm_summary_fixtures.sh
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${REPO_DIR}/tests/data/summary/argnorm"
IMG=$(grep -o "'quay.io/biocontainers/argnorm:[^']*'" "${REPO_DIR}/modules/local/argnorm/main.nf" | tr -d "'")
TMP=$(mktemp -d)
trap 'rm -rf "${TMP}"' EXIT
mkdir -p "${OUT}/args_oap" "${OUT}/deeparg" "${OUT}/card"

# ARGs-OAP 4.0.0 blastout.filtered.txt of sampleA
hit() {  # hit <read> <sseqid> <subtype> <type>
    printf '%s\t%s\t98.0\t50\t50\t300\t1e-20\t90.0\t1.0\tsampleA_tmp_reads\t0.16\t%s\t%s\t%s\t1.0\t0.16\t3.3\t600.0\n' \
        "$1" "$2" "$2" "$3" "$4"
}
{
    printf 'qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore\tqcov\tsample\tscov\tgene\tsubtype\ttype\tcount\tcopy\trpk\trpksum\n'
    n=0
    for spec in "AAG13410 beta_lactam__KPC-2 beta_lactam 3" \
                "HQ451074.1.gene4.p01 beta_lactam__TEM-1 beta_lactam 4" \
                "AAX14397 beta_lactam__mecA beta_lactam 3" \
                "AAX35723 beta_lactam__mecA beta_lactam 1" \
                "AAX19284 vancomycin__vanA vancomycin 2" \
                "ZP_04612557 polymyxin__rosB polymyxin 3" \
                "gb|APB03214.1|ARO:3003980|tetA(48) tetracycline__tetA(48) tetracycline 3"; do
        read -r id sub typ k <<< "${spec}"
        for _ in $(seq 1 "${k}"); do
            n=$((n + 1))
            hit "sampleA_read${n}" "${id}" "${sub}" "${typ}"
        done
    done
} > "${TMP}/sampleA.args_oap.in.tsv"
tail -n +2 "${TMP}/sampleA.args_oap.in.tsv" > "${TMP}/sampleA.args_oap.noheader.tsv"
for s in sampleA sampleB; do
    cp "${REPO_DIR}/tests/data/summary/stages/assembly/05_arg_prediction/contig_level/deeparg/${s}/${s}_contigs_deep_arg.out.mapping.ARG" \
        "${TMP}/${s}.deeparg.in.tsv"
done
docker run --rm --network none -u "$(id -u):$(id -g)" -v "${TMP}:/w" -w /w "${IMG}" sh -c '
    argnorm argsoap -i sampleA.args_oap.noheader.tsv -o sampleA.args_oap.argnorm.tsv
    for s in sampleA sampleB; do argnorm deeparg -i $s.deeparg.in.tsv -o $s.deeparg.argnorm.tsv; done'
install -m 644 "${TMP}/sampleA.args_oap.argnorm.tsv" "${OUT}/args_oap/"
install -m 644 "${TMP}"/sample*.deeparg.argnorm.tsv "${OUT}/deeparg/"

# card/: a hand-written CARD-metadata table (the four columns of the K-MARVEL database's metadata.tsv the summary
# reads: ARO, name, AMR gene family, resistance mechanism) for the AROs the fixtures use, and snps.txt's header. No
# CARD file is copied (roadmap decision RV8: the repository carries no CARD data); tetA(58) (ARO:3003980) is left
# out on purpose (an ARO the staged CARD release lacks)
printf '%s\t%s\t%s\t%s\n' "ARO Accession" "ARO Name" "AMR Gene Family" "Resistance Mechanism" > "${OUT}/card/metadata.tsv"
while IFS='|' read -r aro name family mech; do
    printf '%s\t%s\t%s\t%s\n' "${aro}" "${name}" "${family}" "${mech}" >> "${OUT}/card/metadata.tsv"
done <<'ROWS'
ARO:3000010|vanA|Van ligase;glycopeptide resistance gene cluster|antibiotic target alteration
ARO:3000194|tet(W)|tetracycline-resistant ribosomal protection protein|antibiotic target protection
ARO:3000501|rpoB2|rifamycin-resistant beta-subunit of RNA polymerase (rpoB)|antibiotic target replacement
ARO:3000617|mecA|methicillin resistant PBP2|antibiotic target replacement
ARO:3000838|arlR|major facilitator superfamily (MFS) antibiotic efflux pump|antibiotic efflux
ARO:3000873|TEM-1|TEM beta-lactamase|antibiotic inactivation
ARO:3000980|TEM-117|TEM beta-lactamase|antibiotic inactivation
ARO:3001878|CTX-M-15|CTX-M beta-lactamase|antibiotic inactivation
ARO:3001967|CTX-M-107|CTX-M beta-lactamase|antibiotic inactivation
ARO:3002312|KPC-2|KPC beta-lactamase|antibiotic inactivation
ARO:3002313|KPC-3|KPC beta-lactamase|antibiotic inactivation
ARO:3002724|QnrB10|quinolone resistance protein (qnr)|antibiotic target protection
ARO:3002835|lnuA|lincosamide nucleotidyltransferase (LNU)|antibiotic inactivation
ARO:3003003|CfxA3|CfxA beta-lactamase|antibiotic inactivation
ARO:3003048|rosA|major facilitator superfamily (MFS) antibiotic efflux pump|antibiotic efflux
ARO:3003049|rosB|major facilitator superfamily (MFS) antibiotic efflux pump|antibiotic efflux
ARO:3003206|lsaE|lsa-type ABC-F protein|antibiotic target protection
ARO:3003577|ugd|pmr phosphoethanolamine transferase|antibiotic target alteration
ARO:3003982|LlmA 23S ribosomal RNA methyltransferase|Llm 23S ribosomal RNA methyltransferase|antibiotic target alteration
ARO:3005023|cfrE|Cfr 23S ribosomal RNA methyltransferase|antibiotic target alteration
ARO:3005098|qacL|small multidrug resistance (SMR) antibiotic efflux pump|antibiotic efflux
ARO:3006225|CepA-44|CepA beta-lactamase|antibiotic inactivation
ROWS
printf 'Accession\tName\tModel Type\tParameter Type\tMutations\tCARD Short Name\n' > "${OUT}/card/snps.txt"
# every ARO the argNorm outputs name has a row, except tetA(58)
missing=$(awk -F'\t' 'NR == FNR { have[$1] = 1; next } FNR == 2 { for (i = 1; i <= NF; i++) if ($i == "ARO") c = i }
                       FNR > 2 && !($c in have) && $c != "ARO:3003980" { print $c }' \
    "${OUT}/card/metadata.tsv" "${OUT}"/args_oap/*.tsv "${OUT}"/deeparg/*.tsv | sort -u)
[ -z "${missing}" ] || { echo "AROs without a row in card/metadata.tsv: ${missing}" >&2; exit 1; }
echo "wrote ${OUT}: $(($(wc -l < "${OUT}/card/metadata.tsv") - 1)) metadata rows"
