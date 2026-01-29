# Phase 3: Custom Container Migration Analysis

## Objective
Migrate R scripts from custom `quay.io/ffuentessantander/r_reports:1.1` container to Python scripts using stable biocontainers.

## Modules Using Custom r_reports:1.1 Container

### Reporting Modules (7 modules)
1. **BIN_QUALITY_REPORT** - `bin/r_scripts_temp/Bin_checkm_general_plot.R`
2. **BIN_TAX_REPORT** - `bin/r_scripts_temp/Bins_tax.R`
3. **ARG_BLOBPLOT** - `bin/r_scripts_temp/ARG_blob_plot.R`
4. **ARG_NORM_REPORT** - `bin/r_scripts_temp/Read_arg_norm.R`
5. **ARG_CONTIG_LEVEL_REPORT** - `bin/r_scripts_temp/Contig_arg_unify.R` or similar
6. **BLOBPLOT** - `bin/r_scripts_temp/Blobplot.R`
7. **BIN_SUMMARY** - `bin/r_scripts_temp/Bin_summary.R`

### Database Formatting Modules (7 modules)
8. **FORMAT_SM_DB** - Shell script wrapper
9. **FORMAT_KRAKEN_DB** - Shell script wrapper
10. **FORMAT_BOWTIE_INDEX** - Shell script wrapper
11. **FORMAT_NT_BLAST_DB** - Shell script wrapper
12. **FORMAT_TAXDUMP_FILES** - Shell script wrapper
13. **FORMAT_CHECKM2_DB** - Shell script wrapper
14. **DOWNLOAD_GTDBTK_DB** - Shell script wrapper

## Migration Strategy

### Priority 1: Reporting Modules (High Impact)
These modules generate plots and reports that are critical for pipeline output.

#### 1. BIN_QUALITY_REPORT
**R Script**: `Bin_checkm_general_plot.R`
**Dependencies**: tidyverse, ggplot2
**Functionality**:
- Reads CheckM2 quality reports from multiple binners (MetaBAT2, SemiBin, COMEBin, MetaWRAP)
- Classifies bins as Low/Mid/High quality based on completeness/contamination
- Generates scatter plot (Completeness vs Contamination)
- Outputs: `Total_bins_quality_plot.png`, 3 CSV files

**Python Migration**:
- Use: pandas, matplotlib/seaborn
- Container: `quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0` (pandas + matplotlib)
- Complexity: **Medium** - straightforward data processing and plotting

#### 2. BIN_TAX_REPORT
**R Script**: `Bins_tax.R`
**Dependencies**: tidyverse, rcartocolor
**Functionality**:
- Reads GTDB-Tk taxonomy reports
- Parses GTDB taxonomy strings (d__Bacteria;p__Proteobacteria;...)
- Generates bar plot of MAG counts by Phylum
- Outputs: `MAGs_tax_plot.png`, `MAGs_tax_summary.csv`

**Python Migration**:
- Use: pandas, matplotlib/seaborn
- Container: Same as above
- Complexity: **Medium** - string parsing + bar plot

#### 3. ARG_BLOBPLOT
**R Script**: `ARG_blob_plot.R`
**Dependencies**: tidyverse, rcartocolor, cowplot
**Functionality**:
- Reads contig taxonomy and ARG prediction data
- Creates blob plot (GC vs Coverage, colored by ARG class)
- Creates bar plot of ARG class percentages
- Combines plots with cowplot
- Outputs: `ARG_blob_plot.png`, `ARG_final_blobplot.RDS`

**Python Migration**:
- Use: pandas, matplotlib
- Container: Same as above
- Complexity: **High** - complex multi-panel plot with custom layout

#### 4. ARG_NORM_REPORT
**R Script**: `Read_arg_norm.R`
**Dependencies**: tidyverse
**Functionality**:
- Reads KARGA, KARGVA, and ARGs-OAP results
- Parses complex pipe-delimited gene indices
- Normalizes ARG abundance (CPM, copies per cell)
- Filters by coverage thresholds
- Outputs: `KARGA_norm.csv`, `KARGVA_norm.csv`

**Python Migration**:
- Use: pandas
- Container: `quay.io/biocontainers/pandas:1.5.2`
- Complexity: **Medium** - data processing, no plotting

#### 5. ARG_CONTIG_LEVEL_REPORT
**R Script**: Likely `Contig_arg_unify.R`
**Complexity**: **Medium** - data aggregation

#### 6. BLOBPLOT
**R Script**: `Blobplot.R`
**Complexity**: **High** - similar to ARG_BLOBPLOT

#### 7. BIN_SUMMARY
**R Script**: `Bin_summary.R`
**Complexity**: **Low** - simple data aggregation

### Priority 2: Database Formatting Modules (Low Impact)
These modules are wrappers around shell scripts for database formatting. They don't actually use R - the container is just used for bash execution.

**Migration Strategy**: 
- Replace with minimal container like `ubuntu:22.04` or `quay.io/biocontainers/coreutils:9.1`
- These are simple shell script wrappers, no code migration needed

## Recommended Approach

### Phase 3A: Migrate Simple Reporting Modules
1. **BIN_SUMMARY** (simplest)
2. **ARG_NORM_REPORT** (no plotting)
3. **ARG_CONTIG_LEVEL_REPORT** (data processing)

### Phase 3B: Migrate Plotting Modules
4. **BIN_QUALITY_REPORT** (scatter plot)
5. **BIN_TAX_REPORT** (bar plot)

### Phase 3C: Migrate Complex Plotting Modules
6. **ARG_BLOBPLOT** (multi-panel)
7. **BLOBPLOT** (multi-panel)

### Phase 3D: Update Database Formatting Modules
8-14. Replace r_reports container with minimal bash container

## Container Recommendations

### For Python Reporting Scripts
**Option 1**: `quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0`
- Contains: pandas, matplotlib, numpy, scipy
- Already used in READS_REPORT module
- Stable and well-maintained

**Option 2**: `jupyter/scipy-notebook:latest` (version pinned)
- Contains: pandas, matplotlib, seaborn, numpy, scipy
- Already used in TAXONOMY_REPORT module
- More comprehensive but larger

### For Database Formatting Scripts
**Option**: `ubuntu:22.04` or `quay.io/biocontainers/coreutils:9.1`
- Minimal container for bash scripts
- No R dependencies needed

## Estimated Effort

- **Phase 3A**: 2-3 scripts, ~4-6 hours
- **Phase 3B**: 2 scripts, ~4-6 hours
- **Phase 3C**: 2 scripts, ~6-8 hours
- **Phase 3D**: 7 modules, ~1-2 hours (container swap only)

**Total**: ~15-22 hours of development + testing

## Decision Point

Given the complexity and time investment, we should discuss with the user:

1. **Full Migration**: Migrate all 7 R reporting scripts to Python (~15-20 hours)
2. **Partial Migration**: Migrate simple scripts (3A, 3B), keep complex ones (~8-12 hours)
3. **Container Update Only**: Update to stable R container from Bioconda instead of custom container (~1 hour)

**Recommendation**: Start with **Option 3** (stable R container) for immediate stability, then gradually migrate scripts as needed.
