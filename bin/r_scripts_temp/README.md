# R Scripts Migration - Temporary Folder

This folder contains the original R scripts extracted from `quay.io/ffuentessantander/r_reports:1.1` container.

The migration is complete: every reporting module now calls a Python script in `bin/`
(except phyloseq conversion, kept in R). The original R scripts are **retained here on
purpose** as the reference for visually verifying that the ported scripts produce
equivalent outputs — do not delete this folder until that verification is done.

## Migration Status

**Status Legend:**
- ⏳ Pending
- 🔄 In Progress
- ✅ Migrated
- 🔴 Keep in R

| Module | Original script | Status | Replacement in `bin/` | Notes |
|--------|-----------------|--------|-----------------------|-------|
| reads_report | Report_unify.R | ✅ | `report_unify.py` | |
| arg_norm_report | Read_arg_norm.R | ✅ | `arg_norm_report.py` | |
| taxonomy_report (kraken2 + sourmash) | Tax_unify_report.R | ✅ | `taxonomy_report.py` | One module handles both profilers |
| blobplot | Blobplot.R | ✅ | `blobplot.py` | |
| arg_blobplot | (TBD) | ✅ | `arg_blobplot.py` (+ `blobplot.py`) | |
| bin_summary | Bin_summary.R | ✅ | `bin_summary.py` | |
| bin_quality_report | (TBD) | ✅ | `bin_quality_report.py` | |
| bin_tax_report | (TBD) | ✅ | `bin_tax_report.py` | |
| arg_contig_level_report | (TBD) | ✅ | `arg_contig_level_report.py` | |
| taxonomy_phyloseq | Tax_kraken_to_phyloseq.R | ✅ | `taxonomy_phyloseq.py` | Builds phyloseq-ready tables in Python |
| phyloseq_converter | Tax_kraken_to_phyloseq.R | 🔴 | `tables_to_phyloseq_simple.R` | RDS creation needs the R phyloseq package |

## Verification workflow

1. Run the pipeline (or the relevant module test) with the Python scripts.
2. Compare the outputs (tables/plots) against those produced by the original R scripts in this folder.
3. Once a script's outputs are verified, it no longer needs its R original; remove this folder when all are verified.
