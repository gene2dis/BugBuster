# Phase 3: Python Migration Progress Report

## Status: IN PROGRESS (3/7 reporting modules completed)

## Completed Migrations ✅

### 1. BIN_SUMMARY
- **Original**: `bin/r_scripts_temp/Bin_summary.R` (41 lines)
- **Migrated**: `bin/bin_summary.py` (140 lines)
- **Module**: `modules/local/bin_summary/main.nf` - Updated
- **Container**: `quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0`
- **Dependencies**: pandas
- **Functionality**: Combines bin quality, taxonomy, and depth data into summary table
- **Complexity**: Low - data processing only

### 2. ARG_NORM_REPORT
- **Original**: `bin/r_scripts_temp/Read_arg_norm.R` (165 lines)
- **Migrated**: `bin/arg_norm_report.py` (230 lines)
- **Module**: `modules/local/arg_norm_report/main.nf` - Updated
- **Container**: Same as above
- **Dependencies**: pandas
- **Functionality**: Processes KARGA/KARGVA/ARGs-OAP results, normalizes ARG abundance
- **Complexity**: Medium - complex data parsing and filtering

### 3. ARG_CONTIG_LEVEL_REPORT
- **Original**: `bin/r_scripts_temp/Contig_arg_unify.R` (70 lines)
- **Migrated**: `bin/arg_contig_level_report.py` (120 lines)
- **Module**: `modules/local/arg_contig_level_report/main.nf` - Updated
- **Container**: Same as above
- **Dependencies**: pandas
- **Functionality**: Combines contig taxonomy (blobtools) with ARG predictions (DeepARG)
- **Complexity**: Low - data joining

## Remaining Migrations 🔄

### 4. BIN_QUALITY_REPORT (High Priority)
- **Original**: `bin/r_scripts_temp/Bin_checkm_general_plot.R` (121 lines)
- **Estimated Lines**: ~180 lines Python
- **Module**: `modules/local/bin_quality_report/main.nf`
- **Container**: Same pandas+matplotlib container
- **Dependencies**: pandas, matplotlib
- **Functionality**: 
  - Scatter plot: Completeness vs Contamination
  - Faceted by Raw MAGs vs Refined MAGs
  - Quality classification (Low/Mid/High)
  - 3 CSV outputs + 1 PNG plot
- **Complexity**: Medium - matplotlib scatter plot with facets
- **Estimated Time**: 1-2 hours

### 5. BIN_TAX_REPORT (High Priority)
- **Original**: `bin/r_scripts_temp/Bins_tax.R` (102 lines)
- **Estimated Lines**: ~140 lines Python
- **Module**: `modules/local/bin_tax_report/main.nf`
- **Container**: Same pandas+matplotlib container
- **Dependencies**: pandas, matplotlib
- **Functionality**:
  - Bar plot: MAG counts by Phylum
  - GTDB taxonomy parsing
  - Color palette generation
  - 1 CSV + 1 PNG plot
- **Complexity**: Medium - matplotlib bar plot with custom colors
- **Estimated Time**: 1-2 hours

### 6. ARG_BLOBPLOT (Complex)
- **Original**: `bin/r_scripts_temp/ARG_blob_plot.R` (187 lines)
- **Estimated Lines**: ~250 lines Python
- **Module**: `modules/local/arg_blobplot/main.nf`
- **Container**: Same pandas+matplotlib container
- **Dependencies**: pandas, matplotlib
- **Functionality**:
  - Bubble plot: GC vs Coverage (size=contig length, color=ARG class)
  - Bar plot: ARG class percentages
  - Multi-panel layout (cowplot equivalent)
  - 1 PNG plot + 1 RDS (will be pickle)
- **Complexity**: High - complex multi-panel plot with bubble chart
- **Estimated Time**: 3-4 hours

### 7. BLOBPLOT (Complex)
- **Original**: `bin/r_scripts_temp/Blobplot.R` (estimated ~150 lines)
- **Estimated Lines**: ~200 lines Python
- **Module**: `modules/local/blobplot/main.nf`
- **Container**: Same pandas+matplotlib container
- **Dependencies**: pandas, matplotlib
- **Functionality**:
  - Similar to ARG_BLOBPLOT but for taxonomy
  - Bubble plot with taxonomy coloring
  - Multi-panel layout
- **Complexity**: High - complex multi-panel plot
- **Estimated Time**: 3-4 hours

## Database Formatting Modules (Simple)

These 7 modules only need container updates (no code changes):
- FORMAT_SM_DB
- FORMAT_KRAKEN_DB
- FORMAT_BOWTIE_INDEX
- FORMAT_NT_BLAST_DB
- FORMAT_TAXDUMP_FILES
- FORMAT_CHECKM2_DB
- DOWNLOAD_GTDBTK_DB

**Action**: Replace `r_reports:1.1` with `ubuntu:22.04`
**Estimated Time**: 30 minutes

## Container Used

All Python scripts use the same stable biocontainer:
```
quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0
```

This container includes:
- Python 3.9+
- pandas 1.5.2
- matplotlib 3.6.2
- numpy 1.23.5
- scipy 1.9.3

## Migration Benefits

✅ **Completed Benefits**:
1. Removed dependency on custom `r_reports:1.1` container for 3 modules
2. Using stable, versioned biocontainers
3. Added version tracking to modules
4. Improved error handling and logging
5. More maintainable Python code

## Next Steps

### Option A: Complete All Plotting Modules (~8-10 hours)
Continue with modules 4-7, migrating all plotting scripts to Python with matplotlib.

### Option B: Partial Completion (~2-4 hours)
Complete modules 4-5 (simpler plots), defer complex multi-panel plots (6-7) to future work.

### Option C: Pause for Testing
Test the 3 completed modules before continuing with remaining migrations.

## Recommendation

I recommend **Option A** - completing all migrations now for consistency. The plotting modules are straightforward matplotlib implementations, and completing them ensures:
- No mixed R/Python reporting
- All modules use stable containers
- Consistent code style across pipeline
- Complete Phase 3 closure

**Estimated remaining time**: 8-10 hours for modules 4-7 + 30 minutes for database modules = **8.5-10.5 hours total**

Would you like me to:
1. Continue with all remaining modules (Option A)?
2. Complete simple plotting modules only (Option B)?
3. Pause for testing (Option C)?
