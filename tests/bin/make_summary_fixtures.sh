#!/bin/bash
#
# Regenerate the results-summary fixtures in tests/data/summary/ (used by
# tests/bin/test_results_summary.sh and tests/modules/results_summary.nf.test).
# Two samples, sampleA and sampleB, matching the K-MARVEL fixtures in
# tests/data/kmarvel/ (sampleA: TEM + CTX-M calls, sampleB: no ARGs):
#   sylph_gtdb_{otu_table,tax_table,read_counts}.tsv
#                        bin/taxonomy_phyloseq.py on the sylph-tax captures
#                        tests/data/taxonomy/sample{A,B}_sylph_gtdb.sylphmpa and
#                        their reports (sampleB has no genome, so it is absent
#                        from the tables, as in a real run)
#   Reads_report.csv     TAXONOMY_REPORT layout; read counts are single reads
#                        (2 x the sylph reports' read pairs)
#   sampleA_sylph_gtdb.tsv
#                        sylph profile columns for the three sampleA genomes,
#                        written here: ANI and coverage from the .sylphmpa,
#                        and an ANI interval below 95% for GCA_047199095.1 so
#                        the uncertain-species flag is exercised
#   metadata.tsv         one subject, two visits
#
# Requires python3 with pandas, matplotlib, seaborn and h5py; no network.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
TAX="${REPO_DIR}/tests/data/taxonomy"
OUT="${REPO_DIR}/tests/data/summary"
WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT
mkdir -p "${OUT}"

python3 "${REPO_DIR}/bin/taxonomy_phyloseq.py" \
    --profiler sylph --db-name gtdb --output-dir "${WORK}" --format tables \
    --plot-levels Species --top-n 5 \
    --input-files "${TAX}/sampleA_sylph_gtdb.sylphmpa" "${TAX}/sampleB_sylph_gtdb.sylphmpa" \
                  "${TAX}/sampleA_sylph_gtdb_report.tsv" "${TAX}/sampleB_sylph_gtdb_report.tsv" > /dev/null
for t in otu_table tax_table read_counts; do
    cp "${WORK}/sylph_gtdb_${t}.tsv" "${OUT}/sylph_gtdb_${t}.tsv"
done

cat > "${OUT}/Reads_report.csv" <<'EOF'
Id,Raw reads,Fastp,Clean reads (contaminants removed),Sylph DB,Unclassified,Classified
sampleA,130000,124000,120000,gtdb,0.0,1.0
sampleB,5000,4400,4000,gtdb,1.0,0.0
EOF

printf 'Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tTrue_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tkmers_reassigned\tContig_name\n' > "${OUT}/sampleA_sylph_gtdb.tsv"
printf 'sampleA_R1_clean.fastq.gz\tgtdb_genomes_reps/database/GCA/047/199/115/GCA_047199115.1_genomic.fna.gz\t39.8005\t40.7752\t98.31\t2.234\tNA-NA\tHIGH\tNA-NA\t2\t2.4\t900/1000\t98.31\t0\tcontig_1\n' >> "${OUT}/sampleA_sylph_gtdb.tsv"
printf 'sampleA_R1_clean.fastq.gz\tgtdb_genomes_reps/database/GCA/047/199/095/GCA_047199095.1_genomic.fna.gz\t37.8619\t38.3799\t95.13\t2.125\t93.80-96.40\t1.05\t0.90-1.20\t2\t2.2\t700/1000\t95.13\t0\tcontig_1\n' >> "${OUT}/sampleA_sylph_gtdb.tsv"
printf 'sampleA_R1_clean.fastq.gz\tgtdb_genomes_reps/database/GCF/000/009/045/GCF_000009045.1_genomic.fna.gz\t22.3376\t20.8449\t99.94\t1.254\tNA-NA\tHIGH\tNA-NA\t1\t1.3\t990/1000\t99.94\t0\tcontig_1\n' >> "${OUT}/sampleA_sylph_gtdb.tsv"

printf 'sample\tsubject\tvisit\nsampleA\tS1\tV1\nsampleB\tS1\tV2\n' > "${OUT}/metadata.tsv"

echo "Fixtures written to ${OUT}"
