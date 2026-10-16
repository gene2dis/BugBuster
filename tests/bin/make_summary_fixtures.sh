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
#   rgi/sampleA_rgi_bwt.allele_mapping_data.txt
#                        RGI 6.0.3 bwt allele table (real 24-column header):
#                        KPC-2 on two reference sequences (100% and 60%
#                        coverage, read counts summed), KPC-3 at 30% (below the
#                        80% cutoff), TEM-1 at 95%, CTX-M-15 at 50% (K-MARVEL
#                        calls CTX-M: RGI below cutoff), a WildCARD poxtA from
#                        Blautia producta at 99% and cfrE at 100%
#   rgi/sampleA_rgi_kmer_61mer_analysis.gene.txt
#                        RGI k-mer pathogen of origin: KPC-2 -> K. pneumoniae,
#                        TEM-1 -> E. coli, poxtA -> B. producta, cfrE none
#   rgi/sampleB_rgi_bwt.allele_mapping_data.txt
#                        bin/rgi_bwt_lowcov_fallback.sh placeholder (no hits)
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

mkdir -p "${OUT}/rgi"
H='Reference Sequence\tARO Term\tARO Accession\tReference Model Type\tReference DB\tReference Allele Source\tResistomes & Variants: Observed in Genome(s)\tResistomes & Variants: Observed in Plasmid(s)\tResistomes & Variants: Observed Pathogen(s)\tCompletely Mapped Reads\tMapped Reads with Flanking Sequence\tAll Mapped Reads\tPercent Coverage\tLength Coverage (bp)\tAverage MAPQ (Completely Mapped Reads)\tMate Pair Linkage\tReference Length\tAMR Gene Family\tDrug Class\tResistance Mechanism\tDepth\tSNPs\tConsensus Sequence DNA\tConsensus Sequence Protein\n'
row() {  # row <ref> <term> <aro> <db> <pathogen> <reads> <cov> <len> <family> <drugs> <mechanism>
    printf '%s\t%s\t%s\tprotein homolog model\t%s\tx\tno data\tno data\t%s\t%s\t0\t%s\t%s\t%s\t60\t\t%s\t%s\t%s\t%s\t5.0\tn/a\tACGT\tM\n' \
        "$1" "$2" "$3" "$4" "$5" "$6" "$6" "$7" "$8" "$8" "$9" "${10}" "${11}"
}
{
    printf "${H}"
    row 'ARO:3002312|ID:1|Name:KPC-2' KPC-2 3002312 CARD 'Klebsiella pneumoniae' 30 100.00 882 'KPC beta-lactamase' 'carbapenem;cephalosporin;penicillin beta-lactam' 'antibiotic inactivation'
    row 'Prevalence_Sequence_ID:1|ID:1|Name:KPC-2' KPC-2 3002312 'Resistomes & Variants' 'Klebsiella pneumoniae' 10 60.00 882 'KPC beta-lactamase' 'carbapenem;cephalosporin;penicillin beta-lactam' 'antibiotic inactivation'
    row 'ARO:3002313|ID:2|Name:KPC-3' KPC-3 3002313 CARD 'Klebsiella pneumoniae' 4 30.00 882 'KPC beta-lactamase' 'carbapenem;cephalosporin;penicillin beta-lactam' 'antibiotic inactivation'
    row 'ARO:3000873|ID:3|Name:TEM-1' TEM-1 3000873 CARD 'Escherichia coli' 30 95.00 861 'TEM beta-lactamase' 'monobactam;cephalosporin;penicillin beta-lactam' 'antibiotic inactivation'
    row 'ARO:3001878|ID:4|Name:CTX-M-15' CTX-M-15 3001878 CARD 'Escherichia coli' 6 50.00 876 'CTX-M beta-lactamase' 'cephalosporin' 'antibiotic inactivation'
    row 'Prevalence_Sequence_ID:5|ID:5|Name:poxtA' poxtA 3004470 'Resistomes & Variants' 'Blautia producta' 20 99.00 1632 'Miscellaneous ABC-F subfamily ATP-binding cassette ribosomal protection proteins' 'oxazolidinone antibiotic;phenicol antibiotic;tetracycline antibiotic' 'antibiotic target protection'
    row 'ARO:3005023|ID:6|Name:cfrE' cfrE 3005023 CARD 'Clostridioides difficile' 10 100.00 1050 'Cfr 23S ribosomal RNA methyltransferase' 'lincosamide antibiotic;oxazolidinone antibiotic;phenicol antibiotic' 'antibiotic target alteration'
} > "${OUT}/rgi/sampleA_rgi_bwt.allele_mapping_data.txt"
{
    printf 'ARO term\tMapped reads with kmer DB hits\tCARD*kmer Prediction\tSingle species (chromosome) reads\tSingle species (chromosome or plasmid) reads\tSingle species (plasmid) reads\tSingle species (no genomic info) reads\tSingle genus (chromosome) reads\tSingle genus (chromosome or plasmid) reads\tSingle genus (plasmid) reads\tSingle genus (no genomic info) reads\tPromiscuous plasmid reads\tUnknown taxonomy (chromosome) reads\tUnknown taxonomy (chromosome or plasmid) reads\tUnknown taxonomy (no genomic info) reads\n'
    printf 'KPC-2\t28\tKlebsiella pneumoniae (chromosome or plasmid): 25; Escherichia coli (plasmid): 2; \t\t\t\t\t\t\t\t\t0\t0\t1\t0\n'
    printf 'TEM-1\t9\tEscherichia coli (chromosome or plasmid): 5; \t\t\t\t\t\t\t\t\t0\t0\t4\t0\n'
    printf 'poxtA\t16\tBlautia producta (chromosome): 15; \t\t\t\t\t\t\t\t\t0\t1\t0\t0\n'
    printf 'cfrE\t0\t\t\t\t\t\t\t\t\t\t0\t0\t0\t0\n'
} > "${OUT}/rgi/sampleA_rgi_kmer_61mer_analysis.gene.txt"
(cd "${WORK}" && bash "${REPO_DIR}/bin/rgi_bwt_lowcov_fallback.sh" sampleB_rgi_bwt > /dev/null)
cp "${WORK}/sampleB_rgi_bwt.allele_mapping_data.txt" "${OUT}/rgi/"

echo "Fixtures written to ${OUT}"
