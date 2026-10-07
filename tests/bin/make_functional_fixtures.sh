#!/bin/bash
#
# Generate the committed fixtures used by tests/bin/test_aggregate_functions.sh
# and tests/modules/aggregate_functions.nf.test (functional annotation T4:
# AGGREGATE_FUNCTIONS over featureCounts + eggNOG-mapper v3 output).
#
# Creates, under tests/data/functional/:
#   - sampleA.emapper.annotations  format-faithful eggNOG-mapper v3.0.0-beta6
#     file: variable-length leading '##' block (ctime / version / argv /
#     applied-filters / confidence legend), the exact 22-column '#query'
#     header, one annotated gene (contig_1_1: two KEGG_ko terms, multi-letter
#     COG 'EG', a CAZy family, PFAMs in the real v3 '<name>_<start>_<end>'
#     form (MockPfam twice - a repeated domain must count once - and
#     Mock_dom_2, a name ending in digits; design doc Q15), '-' placeholders, 13-char positional
#     confidence code) with contig_1_2 deliberately absent (unannotated), and
#     the 3 trailing '##' summary lines
#   - sampleA.featureCounts.txt    counts 30/10 over the two genes
#   - sampleB.featureCounts.txt    counts 0/50 (keeps a count-0 gene)
#   - sampleZ.featureCounts.txt    all-zero counts (zero-total-TPM sample)
#   - sampleA.ags.tsv              MICROBECENSUS ags table (T5: GE = 2.0 gives
#     hand-checkable CPGE — RPK 100 and 33.33 over the 300 bp genes -> cpge
#     50 and 16.666667). sampleB/sampleZ deliberately have NO ags fixture:
#     they exercise the TPM-only fallback ('unavailable' in ags_and_ge.tsv)
#   - sampleA.overview.tsv         run_dbcan v5 overview (T6), 8-column
#     layout verified on real 5.2.9 output (incl. the Substrate column the
#     tool's own OVERVIEW_COLUMNS constant omits): contig_1_1 carries calls
#     that exercise the term parsing — '(start-end)' range stripping, '+',
#     ';' and '|' separators (Recommend Results joins with '|' on real
#     output) — with the 'recommended' set {GH5_4, CBM6} deliberately
#     different from the 'any' union {GH5, GH5_4, CBM6} and from eggNOG's
#     CAZy call (GT2); contig_1_2 absent (no CAZyme call)
#   - eggnog_versions.yml          versions.yml shape from EGGNOG_MAPPER_ANNOTATE
#   - dbcan_versions.yml           versions.yml shape from RUN_DBCAN
#
# Gene ids reuse tests/data/gff/test_contig1_genes.gff.gz (contig_1_1 at
# 1-300 +, contig_1_2 at 601-900 -, both 300 bp); the test harnesses copy that
# GFF to per-sample/coassembly names at runtime.
#
# The v3 layout here mirrors the verbatim header of the real 2026-08-29
# acceptance output (design doc Section 4.2). Do not "simplify" the comment
# block: its variable length is part of what the parser must tolerate.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUT_DIR="${SCRIPT_DIR}/../data/functional"
mkdir -p "${OUT_DIR}"
OUT_DIR="$( cd "${OUT_DIR}" && pwd )"

cat > "${OUT_DIR}/sampleA.emapper.annotations" <<'EOF'
## Sat Aug 29 07:21:06 2026
## emapper-3.0.0-beta6
## /usr/local/bin/emapper.py -m no_search --annotate_hits_table sampleA.emapper.seed_orthologs --cpu 8 --data_dir data --output sampleA
##
## applied filters:
##   annot_evalue=0.001
##   annot_score=null
##   tax_scope=auto
## annotation_confidence: one char per annotation field
## confidence codes: h=high m=medium l=low -=not annotated
## confidence field order: Preferred_name GOs EC KEGG_ko KEGG_Pathway KEGG_Module KEGG_Reaction KEGG_rclass BRITE KEGG_TC CAZy BiGG_Reaction PFAMs
#query	seed_ortholog	evalue	score	eggNOG_OGs	tax_ceiling	farthest_donor_lineage	COG_category	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs	annotation_confidence
contig_1_1	1234567.ABC123	2.5e-50	200.0	COG0001@1|root,COG0001@2|Bacteria	2|Bacteria	1|root	EG	mockA	-	-	K00001,K00002	-	-	-	-	-	-	GT2	-	MockPfam_10_95,MockPfam_150_290,Mock_dom_2_5_60	h--h------h-h
## 1 queries scanned
## Total time (seconds): 1.0
## Rate: 1.00 q/s
EOF

cat > "${OUT_DIR}/sampleA.featureCounts.txt" <<'EOF'
# Program:featureCounts v2.1.1; Command:"featureCounts" "-F" "SAF" "-a" "sampleA_genes.saf" "-T" "2" "-p" "--primary" "-o" "sampleA.featureCounts.txt" "sampleA_all_reads.bam"
Geneid	Chr	Start	End	Strand	Length	sampleA_all_reads.bam
contig_1_1	contig_1	1	300	+	300	30
contig_1_2	contig_1	601	900	-	300	10
EOF

cat > "${OUT_DIR}/sampleB.featureCounts.txt" <<'EOF'
# Program:featureCounts v2.1.1; Command:"featureCounts" "-F" "SAF" "-a" "sampleB_genes.saf" "-T" "2" "-p" "--primary" "-o" "sampleB.featureCounts.txt" "sampleB_all_reads.bam"
Geneid	Chr	Start	End	Strand	Length	sampleB_all_reads.bam
contig_1_1	contig_1	1	300	+	300	0
contig_1_2	contig_1	601	900	-	300	50
EOF

cat > "${OUT_DIR}/sampleZ.featureCounts.txt" <<'EOF'
# Program:featureCounts v2.1.1; Command:"featureCounts" "-F" "SAF" "-a" "sampleZ_genes.saf" "-T" "2" "-p" "--primary" "-o" "sampleZ.featureCounts.txt" "sampleZ_all_reads.bam"
Geneid	Chr	Start	End	Strand	Length	sampleZ_all_reads.bam
contig_1_1	contig_1	1	300	+	300	0
contig_1_2	contig_1	601	900	-	300	0
EOF

cat > "${OUT_DIR}/sampleA.ags.tsv" <<'EOF'
sample_id	average_genome_size_bp	genome_equivalents	total_bases
sampleA	3000000	2.0	6000000
EOF

cat > "${OUT_DIR}/eggnog_versions.yml" <<'EOF'
"FUNCTIONAL_ANNOTATION:EGGNOG_MAPPER_ANNOTATE":
    eggnog-mapper: 3.0.0-beta6
    eggnog_db: 7.0.0
EOF

cat > "${OUT_DIR}/sampleA.overview.tsv" <<'EOF'
Gene ID	EC#	dbCAN_hmm	dbCAN_sub	DIAMOND	#ofTools	Recommend Results	Substrate
contig_1_1	3.2.1.4:2|-	GH5(1-95)+CBM6(100-140)	GH5_4(1-95)	GH5;CBM6	3	GH5_4|CBM6	xylan
EOF

cat > "${OUT_DIR}/dbcan_versions.yml" <<'EOF'
"FUNCTIONAL_ANNOTATION:RUN_DBCAN":
    run_dbcan: 5.2.9
    dbcan_db: db_v5-2-9_5-5-2026
EOF

# Per-name copies for tests/modules/aggregate_functions.nf.test: the script
# derives sample ids from filenames, so the GFF must be <id>.gff.gz and the
# coassembly shape needs 'coassembly'-named inputs
cp "${SCRIPT_DIR}/../data/gff/test_contig1_genes.gff.gz" "${OUT_DIR}/sampleA.gff.gz"
cp "${SCRIPT_DIR}/../data/gff/test_contig1_genes.gff.gz" "${OUT_DIR}/coassembly.gff.gz"
sed 's/sampleA/coassembly/g' "${OUT_DIR}/sampleA.emapper.annotations" \
    > "${OUT_DIR}/coassembly.emapper.annotations"
cp "${OUT_DIR}/sampleA.overview.tsv" "${OUT_DIR}/coassembly.overview.tsv"

# Zero-sequence protein FASTA (empty-contig sample shape) for the RUN_DBCAN
# empty-gene-set short-circuit test (-n: no timestamp, deterministic bytes)
gzip -n < /dev/null > "${OUT_DIR}/empty.faa.gz"

echo "Fixtures written under ${OUT_DIR}"
