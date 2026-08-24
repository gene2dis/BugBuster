#!/bin/bash
#
# rgi_bwt_lowcov_fallback.sh <prefix>
#
# Creates placeholder outputs for a sample where `rgi bwt` completed but
# produced no mapping results (too few / no AMR-mapping reads), so RGI_BWT
# still satisfies every non-optional output it declares (audit #8). The
# sorted BAM is deliberately NOT created: there are no alignments, the bam
# emit is optional, and downstream RGI_KMER (which would crash on a fake
# BAM) is simply skipped for the sample.
#
# The allele/gene tables carry the real RGI 6.x headers plus a single NA
# row so RGI_REPORT can parse them alongside real samples.
#
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: rgi_bwt_lowcov_fallback.sh <output_prefix>" >&2
    exit 2
fi

prefix="$1"

allele_cols=(
    "ORF_ID" "Contig" "Start" "Stop" "Orientation" "Cut_Off" "Pass_Bitscore"
    "Best_Hit_ARO" "Best_Identities" "AROMatch" "SNPs_In_Best_Hit_ORTH"
    "Other_Hits" "Unique_Identifier" "Best_Hit_ARO_category" "Best_Resistomes"
    "AROs" "ARO_category" "Resistomes" "Predicted_DNA" "Predicted_Protein"
    "CARD_Protein_Sequence" "Percentage_Length_of_CARD_Protein" "ID"
    "Model_ID" "Nudged" "Note" "Other_Hit_Accession"
)
gene_cols=(
    "ORF_ID" "ARO Term" "ARO Accession" "Reference Model Type" "Reference DB"
    "Alleles with Mapped Reads" "Reference Allele"
    "%Coverage of Reference Allele" "Minimum Bidirectional Coverage"
    "Average Bidirectional Coverage" "%Identity to Reference Allele"
    "Antibiotic" "Class" "Resistance Mechanism" "AMR Gene Family" "Drug Class"
)

write_na_table() {
    # write_na_table <outfile> <col>...
    local outfile="$1"
    shift
    local header na
    header=$(printf '%s\t' "$@")
    na=$(printf 'NA\t%.0s' "$@")
    printf '%s\n%s\n' "${header%$'\t'}" "${na%$'\t'}" > "$outfile"
}

write_na_table "${prefix}.allele_mapping_data.txt" "${allele_cols[@]}"
write_na_table "${prefix}.gene_mapping_data.txt" "${gene_cols[@]}"

echo "No mapping statistics available - insufficient reads mapped to CARD" \
    > "${prefix}.overall_mapping_stats.txt"
echo "No mapping statistics available - insufficient reads mapped to CARD" \
    > "${prefix}.reference_mapping_stats.txt"
