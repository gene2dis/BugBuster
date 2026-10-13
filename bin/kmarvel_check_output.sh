#!/bin/bash
#
# kmarvel_check_output.sh <kmarvel-stdout-log> <kmarvel-tsv>
#
# Turns two K-MARVEL 1.0.0 (commit 7108d2d) exit-0 outcomes into what the
# KMARVEL module needs:
#   - zero ARGs found: kmarvel writes its TSV with no rows and no header (the
#     csv writer emits the header with the first record). Writes the
#     header-only table, so KMARVEL_REPORT keeps the sample as all zeros.
#   - protein BLAST failed: kmarvel prints "Protein BLAST failed ..." and
#     still exits 0, writing unfiltered results without the blast_* columns.
#     With ARGs found that is a broken result: exit 1 so the task fails.
#     (With zero ARGs blastp runs on an empty query; nothing is lost.)
# A log without the "Found N Antimicrobial Resistance Genes" line means
# kmarvel stopped early or its output changed: exit 1.
#
# Remove with the re-pin once upstream exits non-zero on a BLAST failure
# and always writes the TSV header.
#
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: kmarvel_check_output.sh <kmarvel_stdout_log> <kmarvel_tsv>" >&2
    exit 2
fi

log="$1"
tsv="$2"

# FlatFinalFilteredData, in field order (src/main.rs at 7108d2d)
cols=(
    aro_full_name aro_id gene_family drug_class resistance_mechanism model
    identity assembly_coverage supporting_reads kmer_ratio median_kmer_depth
    global_alignment_score local_alignment_score composite_score
    blast_predicted_aro_name blast_predicted_aro_id blast_predicted_gene_family
    blast_predicted_drug_class blast_predicted_resistance_mechanism
    blast_e_value blast_bit_score blast_assembled_sequence
    blast_reference_sequence blast_cigar_string blast_mutation_data
    assembled_sequence reference_sequence cigar_string mutation_data
    start end read_data
)

found=$(sed -n 's/^\[INFO\] Found \([0-9][0-9]*\) Antimicrobial Resistance Genes (ARGs) for .*/\1/p' "${log}" | tail -1)
if [ -z "${found}" ]; then
    echo "ERROR: no 'Found N Antimicrobial Resistance Genes' line in ${log}: K-MARVEL did not finish" >&2
    exit 1
fi

if [ "${found}" -eq 0 ]; then
    echo "K-MARVEL found no ARGs: writing a header-only ${tsv}"
    (IFS=$'\t'; printf '%s\n' "${cols[*]}") > "${tsv}"
    exit 0
fi

if grep -q '^Protein BLAST \(failed\|terminated\)' "${log}"; then
    echo "ERROR: K-MARVEL found ${found} ARGs but its protein BLAST step failed (see ${log}); the results lack the BLAST filter" >&2
    exit 1
fi

if [ ! -s "${tsv}" ] || ! head -1 "${tsv}" | tr '\t' '\n' | grep -qx 'blast_e_value'; then
    echo "ERROR: ${tsv} is missing or has no blast_e_value column although K-MARVEL found ${found} ARGs" >&2
    exit 1
fi
