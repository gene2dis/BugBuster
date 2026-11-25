#!/bin/bash
#
# Tests for bin/tables_to_phyloseq_simple.R, PHYLOSEQ_CONVERTER's script
# (--create_phyloseq_rds; roadmap Phase 25, RV25; review findings R6-15,
# R8-02), run in the image the module pins and read back with phyloseq:
#   - the RDS holds every taxon and sample of the TAXONOMY_PHYLOSEQ tables
#   - taxon names with an apostrophe or a '#', and a taxon or sample called
#     NA, are kept as they are (roadmap Phase 20, RV6; review R6-02: the
#     read.table defaults merged or cut such rows, or made NA missing)
#   - a taxonomy table whose rows do not match the OTU table stops the script
#     with its message, and so do too few arguments
#
# Requires docker.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"

# The image the module pins, so the test cannot drift from it
PHYLOSEQ_IMG=$(grep -o "container 'quay.io/biocontainers/bioconductor-phyloseq:[^']*'" \
    "${REPO_DIR}/modules/local/phyloseq_converter/main.nf" | sed "s/^container '//; s/'\$//")
if [ -z "${PHYLOSEQ_IMG}" ]; then
    echo "ERROR: could not read the phyloseq image from modules/local/phyloseq_converter/main.nf" >&2
    exit 1
fi

TMP_DIR=$(mktemp -d)
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

PASS=0
FAIL=0
ok() { echo "✓ $1"; PASS=$((PASS + 1)); }
ko() { echo "✗ $1"; FAIL=$((FAIL + 1)); }

in_image() {
    # in_image <command...>: run in the pinned image with bin/ and the work dir
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" -v "${TMP_DIR}:/w" -w /w \
        "${PHYLOSEQ_IMG}" "$@"
}

echo "=== tables_to_phyloseq_simple.R test suite ==="
echo "Image: ${PHYLOSEQ_IMG}"

# Tables as taxonomy_phyloseq.py writes them (pandas to_csv, tab-separated,
# index in the first column, a field quoted only when it holds a double quote)
printf '\tsample1\tNA\n562\t610.0\t380.0\nNA\t5.0\t0.0\n1280\t0.0\t420.0\n' > "${TMP_DIR}/otu.tsv"
printf '\tKingdom\tGenus\tSpecies\n562\tBacteria\tEscherichia\tEscherichia coli\nNA\tBacteria\tCandidatus O'"'"'Brienella\tBacterium #2\n1280\tBacteria\tStaphylococcus\t"Staphylococcus ""aureus"""\n' > "${TMP_DIR}/tax.tsv"
printf '\tsample_name\nsample1\tsample1\nNA\tNA\n' > "${TMP_DIR}/meta.tsv"

if in_image Rscript /pipeline_bin/tables_to_phyloseq_simple.R otu.tsv tax.tsv meta.tsv out.RDS > "${TMP_DIR}/run.log" 2>&1; then
    ok "builds the RDS from three tables"
else
    ko "builds the RDS from three tables: $(tail -3 "${TMP_DIR}/run.log")"
fi

# Read the object back: counts, the names as written, and one abundance
cat > "${TMP_DIR}/read.R" <<'EOF'
suppressPackageStartupMessages(library(phyloseq))
ps <- readRDS("out.RDS")
cat("ntaxa", ntaxa(ps), "\n")
cat("nsamples", nsamples(ps), "\n")
cat("taxa", paste(taxa_names(ps), collapse = "|"), "\n")
cat("samples", paste(sample_names(ps), collapse = "|"), "\n")
tt <- as(tax_table(ps), "matrix")
cat("genus_NA", tt["NA", "Genus"], "\n")
cat("species_NA", tt["NA", "Species"], "\n")
cat("species_1280", tt["1280", "Species"], "\n")
cat("count_1280_NA", as(otu_table(ps), "matrix")["1280", "NA"], "\n")
EOF
in_image Rscript read.R > "${TMP_DIR}/read.log" 2>&1 || true
field() { sed -n "s/^$1 //p" "${TMP_DIR}/read.log" | sed 's/ $//'; }
[ "$(field ntaxa)" = "3" ] && ok "three taxa" || ko "three taxa: $(field ntaxa) ($(cat "${TMP_DIR}/read.log"))"
[ "$(field nsamples)" = "2" ] && ok "two samples" || ko "two samples: $(field nsamples)"
[ "$(field taxa)" = "562|NA|1280" ] && ok "a taxon called NA keeps its name" || ko "taxa: $(field taxa)"
[ "$(field samples)" = "sample1|NA" ] && ok "a sample called NA keeps its name" || ko "samples: $(field samples)"
[ "$(field genus_NA)" = "Candidatus O'Brienella" ] && ok "a name with an apostrophe is read whole" || ko "genus with apostrophe: $(field genus_NA)"
[ "$(field species_NA)" = "Bacterium #2" ] && ok "a name with '#' is not cut" || ko "species with '#': $(field species_NA)"
[ "$(field species_1280)" = 'Staphylococcus "aureus"' ] && ok "a quoted field with doubled quotes is unquoted" || ko "quoted species: $(field species_1280)"
[ "$(field count_1280_NA)" = "420" ] && ok "abundances stay with their taxon and sample" || ko "count 1280/NA: $(field count_1280_NA)"

# A taxonomy table that does not describe the same taxa
printf '\tKingdom\tGenus\tSpecies\n562\tBacteria\tEscherichia\tEscherichia coli\n' > "${TMP_DIR}/tax_short.tsv"
if in_image Rscript /pipeline_bin/tables_to_phyloseq_simple.R otu.tsv tax_short.tsv meta.tsv bad.RDS > "${TMP_DIR}/bad.log" 2>&1; then
    ko "row-count mismatch stops the script: it succeeded"
elif grep -q "the taxonomy table has 1 rows and the OTU table 3" "${TMP_DIR}/bad.log" && [ ! -e "${TMP_DIR}/bad.RDS" ]; then
    ok "row-count mismatch stops the script with its message, no RDS"
else
    ko "row-count mismatch: wrong message ($(tail -2 "${TMP_DIR}/bad.log"))"
fi

if in_image Rscript /pipeline_bin/tables_to_phyloseq_simple.R otu.tsv > "${TMP_DIR}/args.log" 2>&1; then
    ko "too few arguments: the script succeeded"
elif grep -q "Error: Insufficient arguments" "${TMP_DIR}/args.log"; then
    ok "too few arguments stop the script with the usage"
else
    ko "too few arguments: wrong message ($(tail -2 "${TMP_DIR}/args.log"))"
fi

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
