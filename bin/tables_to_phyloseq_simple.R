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

# Read OTU table
otu_data <- read.table(
  otu_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Read taxonomy table
tax_data <- read.table(
  tax_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Read sample metadata
sample_data <- read.table(
  sample_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

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
