# BugBuster Pipeline Refactoring Plan

## Overview
This document outlines the detailed plan to fix critical issues in the BugBuster pipeline after the assembly and binning modules refactoring.

## Issues Identified

### 1. Output Structure Misalignment
**Problem**: The pipeline outputs don't follow the structure defined in `config/modules.config`

**Current Issues**:
- Assembly outputs should go to `03_assembly/per_sample/{sample_id}` or `03_assembly/coassembly`
- Binning outputs should go to `04_binning/per_sample/{sample_id}` or `04_binning/coassembly`
- Some modules may not be respecting the publishDir configurations

**Files to Check**:
- `config/modules.config` (lines 224-460) - publishDir definitions
- `subworkflows/local/assembly.nf` - MEGAHIT, BBMAP outputs
- `subworkflows/local/binning.nf` - METABAT2, SEMIBIN, COMEBIN, METAWRAP outputs

### 2. Custom Container Dependencies
**Problem**: Multiple modules use custom container `quay.io/ffuentessantander/r_reports:1.1` with R scripts embedded in the container

**Affected Modules** (9 total):
1. `BIN_QUALITY_REPORT` - uses `/mnt/Bin_checkm_general_plot.R`
2. `BIN_TAX_REPORT` - uses `/mnt/Bins_tax.R`
3. `ARG_BLOBPLOT` - uses R script from container
4. `ARG_NORM_REPORT` - uses R script from container
5. `BLOBPLOT` - uses R script from container
6. `ARG_CONTIG_LEVEL_REPORT` - uses R script from container
7. `BIN_SUMMARY` - uses R script from container
8. `FORMAT_SM_DB` - uses bash script from container
9. Database formatting processes (FORMAT_KRAKEN_DB, FORMAT_BOWTIE_INDEX, etc.)

**R Scripts Available in bin/r_scripts_temp/**:
- ARG_blob_plot.R
- Bin_checkm_general_plot.R
- Bin_checkm_per_sample_plot.R
- Bin_summary.R
- Bins_tax.R
- Blobplot.R
- Contig_arg_unify.R
- Read_arg_norm.R
- Tax_kraken_to_phyloseq.R
- Tax_sourmash_to_phyloseq.R
- Tax_unify_report.R

**Python Scripts Already in bin/**:
- report_unify.py (used by READS_REPORT with stable container)
- taxonomy_report.py (used by TAXONOMY_REPORT with stable container)
- taxonomy_phyloseq.py (used by TAXONOMY_PHYLOSEQ)

**Migration Strategy**:
- Option A: Migrate R scripts to Python (preferred for easier dependency management)
- Option B: Use stable R container (e.g., `rocker/tidyverse`) and move scripts to bin/

### 3. Performance Inefficiencies
**Problem**: Some steps are taking too long to run

**Areas to Investigate**:
1. **CPU Allocation**:
   - MEGAHIT: `process_medium` (8 CPUs, 36GB RAM)
   - METABAT2: `process_medium` (8 CPUs, 36GB RAM)
   - SEMIBIN: `process_medium` (8 CPUs, 36GB RAM)
   - COMEBIN: `process_high` (16 CPUs, 72GB RAM)
   - BOWTIE2: `process_medium` (8 CPUs, 36GB RAM)

2. **Parallelization**:
   - Check if multiple samples can run in parallel
   - Review executor settings in `conf/base.config`
   - Consider maxForks settings for resource-intensive processes

3. **Known Issues from TODO.md**:
   - Mapping steps taking too long
   - Need to optimize annotation modules (deeparg, prodigal, metaceberus)
   - Other modules (autometa, vamb) may have issues

### 4. Incorrect FASTQ File Usage (CRITICAL)
**Problem**: Downstream processes use host-cleaned reads instead of PhiX-cleaned reads

**Current Flow in `subworkflows/local/qc.nf`**:
```
FASTP → QFILTER → BOWTIE2_HOST (remove host) → BOWTIE2_PHIX (remove PhiX)
                                      ↓                           ↓
                              ch_host_clean.reads         ch_phix_clean.reads
```

**Lines 104-105 in qc.nf**:
```groovy
ch_clean_reads = ch_host_clean.reads              // ❌ WRONG - uses host-cleaned only
ch_clean_reads_coassembly = ch_host_clean.reads_coassembly  // ❌ WRONG
```

**Should be**:
```groovy
ch_clean_reads = ch_phix_clean.reads              // ✓ CORRECT - uses PhiX-cleaned
ch_clean_reads_coassembly = ch_phix_clean.reads_coassembly  // ✓ CORRECT
```

**Impact**: All downstream processes receive reads that still contain PhiX contamination:
- Assembly (MEGAHIT)
- Binning (METABAT2, SEMIBIN, COMEBIN)
- Taxonomic profiling (KRAKEN2, SOURMASH)
- ARG prediction (KARGA, KARGVA, ARGS_OAP)

## Detailed Action Plan

### Phase 1: Critical Fix - FASTQ File Routing (HIGHEST PRIORITY)
**Estimated Time**: 5 minutes
**Files**: `subworkflows/local/qc.nf`

1. Fix lines 104-105 to use `ch_phix_clean` instead of `ch_host_clean`
2. Verify the change propagates correctly to all downstream processes

### Phase 2: Output Structure Validation and Fixes
**Estimated Time**: 30-45 minutes
**Files**: Multiple module files and config

1. Audit all assembly and binning module publishDir settings
2. Compare against `config/modules.config` specifications
3. Fix any misalignments in:
   - `modules/local/megahit/main.nf`
   - `modules/local/bbmap/main.nf`
   - `modules/local/metabat2/main.nf`
   - `modules/local/semibin/main.nf`
   - `modules/local/comebin/main.nf`
   - `modules/local/metawrap/main.nf`
   - `modules/local/checkm2/main.nf`
   - `modules/local/gtdb-tk/main.nf`

### Phase 3: Container Migration
**Estimated Time**: 2-3 hours
**Files**: 9 module files + new Python scripts

**Approach**: Migrate R scripts to Python for better reproducibility

1. **Priority 1 - Binning Reports** (most critical):
   - Migrate `Bin_checkm_general_plot.R` → `bin/bin_quality_report.py`
   - Migrate `Bins_tax.R` → `bin/bin_tax_report.py`
   - Update `modules/local/bin_quality_report/main.nf`
   - Update `modules/local/bin_tax_report/main.nf`
   - Use stable container: `quay.io/biocontainers/mulled-v2-*` with pandas, matplotlib, seaborn

2. **Priority 2 - ARG Reports**:
   - Migrate `Read_arg_norm.R` → `bin/arg_norm_report.py`
   - Migrate `ARG_blob_plot.R` → `bin/arg_blobplot.py`
   - Migrate `Contig_arg_unify.R` → `bin/arg_contig_level_report.py`
   - Update corresponding modules

3. **Priority 3 - Taxonomy/Blob Reports**:
   - Migrate `Blobplot.R` → `bin/blobplot.py`
   - Update `modules/local/blobplot/main.nf`

4. **Priority 4 - Bin Summary**:
   - Migrate `Bin_summary.R` → `bin/bin_summary.py`
   - Update `modules/local/bin_summary/main.nf`

5. **Database Formatting** (if needed):
   - Review if custom scripts are truly needed
   - Consider using standard biocontainers tools

### Phase 4: Performance Optimization
**Estimated Time**: 1-2 hours
**Files**: `nextflow.config`, `conf/base.config`, module files

1. **Review and Adjust CPU Allocations**:
   - Analyze typical resource usage patterns
   - Consider increasing MEGAHIT to `process_high` if assemblies are large
   - Review BOWTIE2_SAMTOOLS parallelization

2. **Optimize Parallelization**:
   - Add `maxForks` directives for memory-intensive processes
   - Review executor settings for better throughput
   - Consider process-specific queue configurations

3. **Specific Optimizations**:
   - BOWTIE2 alignment: Review if streaming to SAMtools is optimal
   - METABAT2: Check if depth calculation can be parallelized
   - Review if any processes can be combined to reduce I/O

4. **Address TODO.md Items**:
   - Investigate slow mapping steps
   - Review annotation module configurations

### Phase 5: Testing and Validation
**Estimated Time**: Variable
**Files**: Test datasets

1. Create test script to validate:
   - Correct FASTQ files are used at each step
   - Output structure matches config
   - All containers are stable/public
   - Performance improvements are measurable

2. Run with test profile:
   ```bash
   nextflow run main.nf -profile test,docker
   ```

3. Validate outputs:
   - Check directory structure
   - Verify file contents
   - Compare performance metrics

## Implementation Order

1. ✅ **CRITICAL**: Fix FASTQ file routing (Phase 1)
2. ✅ Validate output structure (Phase 2)
3. ✅ Migrate custom containers (Phase 3)
4. ✅ Optimize performance (Phase 4)
5. ✅ Test and validate (Phase 5)

## Success Criteria

- [ ] All downstream processes use PhiX-cleaned reads
- [ ] Output directory structure matches `config/modules.config`
- [ ] No custom/unstable containers in use
- [ ] All scripts are in `bin/` directory
- [ ] Performance improvements documented
- [ ] Test profile runs successfully
- [ ] Documentation updated

## Notes

- Keep original R scripts in `bin/r_scripts_temp/` as reference
- Document any breaking changes
- Update README with new requirements
- Consider adding integration tests
