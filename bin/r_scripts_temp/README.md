# R Scripts Migration - Temporary Folder

This folder contains the original R scripts extracted from `quay.io/ffuentessantander/r_reports:1.1` container.

## Scripts to Migrate

### Reporting/Analysis Scripts
- `Report_unify.R` → used by `reads_report` module
- `Read_arg_norm.R` → used by `arg_norm_report` module
- `Tax_unify_report.R` → used by `tax_report_kraken2` and `tax_report_sourmash` modules
- `Blobplot.R` → used by `blobplot` module
- `Bin_summary.R` → used by `bin_summary` module
- `Tax_kraken_to_phyloseq.R` → used by `kraken_to_phyloseq` module
- Additional scripts for other modules (arg_blobplot, bin_quality_report, bin_tax_report, etc.)

## Migration Status

Scripts will be migrated incrementally, module by module.

**Status Legend:**
- ⏳ Pending
- 🔄 In Progress
- ✅ Migrated
- 🔴 Keep in R

| Module | Script | Status | Target Language | Notes |
|--------|--------|--------|-----------------|-------|
| reads_report | Report_unify.R | ⏳ | TBD | |
| arg_norm_report | Read_arg_norm.R | ⏳ | TBD | |
| tax_report_kraken2 | Tax_unify_report.R | ⏳ | TBD | |
| tax_report_sourmash | Tax_unify_report.R | ⏳ | TBD | |
| blobplot | Blobplot.R | ⏳ | TBD | |
| arg_blobplot | (TBD) | ⏳ | TBD | |
| bin_summary | Bin_summary.R | ⏳ | TBD | |
| bin_quality_report | (TBD) | ⏳ | TBD | |
| bin_tax_report | (TBD) | ⏳ | TBD | |
| kraken_to_phyloseq | Tax_kraken_to_phyloseq.R | ⏳ | TBD | Likely keep R (phyloseq) |
| sourmash_to_phyloseq | (TBD) | ⏳ | TBD | Likely keep R (phyloseq) |
| arg_contig_level_report | (TBD) | ⏳ | TBD | |

## Workflow

1. Add original R scripts to this folder
2. Analyze script for migration (per user request)
3. Create migrated version in `bin/` (Python or R)
4. Update module to use new script and container
5. Test and validate
6. Mark as complete in table above
