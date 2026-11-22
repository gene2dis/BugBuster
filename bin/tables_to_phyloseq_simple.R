#!/usr/bin/env Rscript
#
# Convert TSV tables to R phyloseq object (simplified version without optparse)
#
# Usage: Rscript tables_to_phyloseq_simple.R <otu_table.tsv> <tax_table.tsv> <sample_metadata.tsv> <output.RDS>
#
# Author: BugBuster Development Team
# License: MIT

suppressPackageStartupMessages({
  library(phyloseq)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Error: Insufficient arguments\n")
  cat("Usage: Rscript tables_to_phyloseq_simple.R <otu_table.tsv> <tax_table.tsv> <sample_metadata.tsv> <output.RDS>\n")
  quit(status = 1)
}

otu_file <- args[1]
tax_file <- args[2]
sample_file <- args[3]
output_file <- args[4]

# Validate input files exist
if (!file.exists(otu_file)) {
  cat("Error: OTU table file not found:", otu_file, "\n")
  quit(status = 1)
}
if (!file.exists(tax_file)) {
  cat("Error: Taxonomy table file not found:", tax_file, "\n")
  quit(status = 1)
}
if (!file.exists(sample_file)) {
  cat("Error: Sample metadata file not found:", sample_file, "\n")
  quit(status = 1)
}

cat("Reading input files...\n")

# The tables are written by taxonomy_phyloseq.py with pandas' to_csv, which
# wraps a field in double quotes only when it holds one (doubled inside), a tab
# or a line break, and never quotes an apostrophe or a '#'. read.table's
# default quote ("'" as well) and comment ("#") characters merged or cut rows
# of taxon names holding an apostrophe or a '#', with no error, and its default
# na.strings made a taxon or sample called NA missing (roadmap Phase 20, RV6;
# review finding R6-02). The reads below take double quotes as pandas writes
# them and every other character as it is.

# Read OTU table
otu_data <- read.table(
  otu_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\"",
  comment.char = "",
  na.strings = character(0)
)

# Read taxonomy table
tax_data <- read.table(
  tax_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\"",
  comment.char = "",
  na.strings = character(0)
)

# Read sample metadata
sample_data <- read.table(
  sample_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\"",
  comment.char = "",
  na.strings = character(0)
)

if (nrow(tax_data) != nrow(otu_data)) {
  cat("Error: the taxonomy table has", nrow(tax_data), "rows and the OTU table", nrow(otu_data),
      "- the tables do not describe the same taxa\n")
  quit(status = 1)
}

cat("Creating phyloseq object...\n")

# Convert to phyloseq components
otu_matrix <- as.matrix(otu_data)
tax_matrix <- as.matrix(tax_data)

# Create phyloseq object
ps <- phyloseq(
  otu_table(otu_matrix, taxa_are_rows = TRUE),
  tax_table(tax_matrix),
  sample_data(sample_data)
)

# Save phyloseq object
cat("Saving phyloseq object to:", output_file, "\n")
saveRDS(ps, file = output_file)

cat("Successfully created phyloseq object with:\n")
cat("  -", ntaxa(ps), "taxa\n")
cat("  -", nsamples(ps), "samples\n")
cat("  -", ncol(tax_table(ps)), "taxonomic ranks\n")
cat("\nDone!\n")
