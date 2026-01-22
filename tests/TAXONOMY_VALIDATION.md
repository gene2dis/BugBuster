# Taxonomy Workflow Validation Guide

This document provides validation procedures for the refactored taxonomy workflow.

## Overview

The refactored taxonomy workflow replaces 5 R-based modules with 2 Python-based core modules (+1 optional R converter), improving performance, maintainability, and interoperability.

## Validation Checklist

### 1. Script Functionality Tests

#### Test taxonomy_report.py
```bash
# Test help
python3 bin/taxonomy_report.py --help

# Test with sample Kraken2 data
python3 bin/taxonomy_report.py \
    --profiler kraken2 \
    --reports tests/data/taxonomy/sample1_kraken_report.tsv \
    --reads-report tests/data/taxonomy/Reads_report.csv \
    --db-name silva \
    --output-dir test_output/

# Expected outputs:
# - test_output/Reads_report.csv (updated)
# - test_output/Kraken_plot.png
```

#### Test taxonomy_phyloseq.py
```bash
# Test help
python3 bin/taxonomy_phyloseq.py --help

# Test with Kraken2 BIOM data
python3 bin/taxonomy_phyloseq.py \
    --profiler kraken2 \
    --input-files tests/data/taxonomy/sample1.biom \
    --db-name silva \
    --output-dir test_output/ \
    --format both \
    --plot-levels Phylum,Family,Genus,Species \
    --top-n 10

# Expected outputs:
# - test_output/kraken2_silva_otu_table.tsv
# - test_output/kraken2_silva_tax_table.tsv
# - test_output/kraken2_silva_sample_metadata.tsv
# - test_output/kraken2_silva_phyloseq_data.h5
# - test_output/plots/*.png (4 plots)
```

#### Test tables_to_phyloseq.R
```bash
# Test help
Rscript bin/tables_to_phyloseq.R --help

# Test conversion
Rscript bin/tables_to_phyloseq.R \
    --otu-table test_output/kraken2_silva_otu_table.tsv \
    --tax-table test_output/kraken2_silva_tax_table.tsv \
    --sample-data test_output/kraken2_silva_sample_metadata.tsv \
    --output test_output/kraken2_silva_phyloseq.RDS \
    --validate

# Expected output:
# - test_output/kraken2_silva_phyloseq.RDS
```

### 2. Module Integration Tests

#### Test TAXONOMY_REPORT module
```bash
# Run basic module test
nextflow run tests/modules/taxonomy/taxonomy_report.nf.test
```

#### Test TAXONOMY_PHYLOSEQ module
```bash
# Run basic module test
nextflow run tests/modules/taxonomy/taxonomy_phyloseq.nf.test
```

#### Test PHYLOSEQ_CONVERTER module
```bash
# Run basic module test
nextflow run tests/modules/taxonomy/phyloseq_converter.nf.test
```

### 3. Subworkflow Integration Tests

#### Test complete Kraken2 pathway
```bash
nextflow run main.nf \
    -profile test,docker \
    --taxonomic_profiler kraken2 \
    --taxonomy_plot_levels "Phylum,Family,Genus" \
    --taxonomy_top_n_taxa 10 \
    --create_phyloseq_rds false
```

**Expected outputs**:
- `results/reports/read_level/Reads_report.csv`
- `results/reports/read_level/taxonomy/Kraken_plot.png`
- `results/reports/read_level/taxonomy/kraken2_*_otu_table.tsv`
- `results/reports/read_level/taxonomy/kraken2_*_tax_table.tsv`
- `results/reports/read_level/taxonomy/kraken2_*_sample_metadata.tsv`
- `results/reports/read_level/taxonomy/kraken2_*_phyloseq_data.h5`
- `results/reports/read_level/taxonomy/plots/*.png`

#### Test complete Sourmash pathway
```bash
nextflow run main.nf \
    -profile test,docker \
    --taxonomic_profiler sourmash \
    --taxonomy_plot_levels "Phylum,Family,Genus,Species" \
    --taxonomy_top_n_taxa 15 \
    --create_phyloseq_rds false
```

**Expected outputs**:
- `results/reports/read_level/Reads_report.csv`
- `results/reports/read_level/taxonomy/sourmash_tax_classified_reads.png`
- `results/reports/read_level/taxonomy/sourmash_*_otu_table.tsv`
- `results/reports/read_level/taxonomy/sourmash_*_tax_table.tsv`
- `results/reports/read_level/taxonomy/sourmash_*_sample_metadata.tsv`
- `results/reports/read_level/taxonomy/sourmash_*_phyloseq_data.h5`
- `results/reports/read_level/taxonomy/plots/*.png`

#### Test with R phyloseq generation
```bash
nextflow run main.nf \
    -profile test,docker \
    --taxonomic_profiler kraken2 \
    --create_phyloseq_rds true
```

**Additional expected output**:
- `results/reports/read_level/taxonomy/kraken2_*_phyloseq.RDS`

### 4. Output Validation

#### Validate TSV table structure
```python
import pandas as pd

# Load OTU table
otu = pd.read_csv('results/.../otu_table.tsv', sep='\t', index_col=0)
assert otu.shape[0] > 0, "OTU table has no taxa"
assert otu.shape[1] > 0, "OTU table has no samples"
assert (otu >= 0).all().all(), "OTU table has negative values"

# Load taxonomy table
tax = pd.read_csv('results/.../tax_table.tsv', sep='\t', index_col=0)
assert tax.shape[0] == otu.shape[0], "OTU and tax tables have different taxa counts"
expected_ranks = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
assert all(rank in tax.columns for rank in expected_ranks), "Missing taxonomic ranks"

# Load sample metadata
meta = pd.read_csv('results/.../sample_metadata.tsv', sep='\t', index_col=0)
assert meta.shape[0] == otu.shape[1], "Metadata and OTU table have different sample counts"
assert 'profiler' in meta.columns, "Missing profiler column"
assert 'database' in meta.columns, "Missing database column"

print("✓ All table validations passed")
```

#### Validate HDF5 structure
```python
import h5py

with h5py.File('results/.../phyloseq_data.h5', 'r') as f:
    assert 'otu_table' in f, "Missing OTU table in HDF5"
    assert 'tax_table' in f, "Missing taxonomy table in HDF5"
    assert 'metadata' in f, "Missing metadata in HDF5"
    
    otu_shape = f['otu_table'].shape
    tax_shape = f['tax_table'].shape
    assert otu_shape[0] == tax_shape[0], "OTU and tax tables have different taxa counts"
    
    print(f"✓ HDF5 validation passed: {otu_shape[0]} taxa × {otu_shape[1]} samples")
```

#### Validate R phyloseq object (if generated)
```R
library(phyloseq)

ps <- readRDS('results/.../phyloseq.RDS')

# Check phyloseq object structure
stopifnot(class(ps) == "phyloseq")
stopifnot(ntaxa(ps) > 0)
stopifnot(nsamples(ps) > 0)

# Check components
stopifnot(!is.null(otu_table(ps)))
stopifnot(!is.null(tax_table(ps)))

# Check taxonomic ranks
ranks <- rank_names(ps)
expected_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
stopifnot(all(expected_ranks %in% ranks))

cat("✓ R phyloseq object validation passed\n")
cat(sprintf("  Taxa: %d\n", ntaxa(ps)))
cat(sprintf("  Samples: %d\n", nsamples(ps)))
```

### 5. Performance Validation

#### Compare execution times
```bash
# Old workflow (if available for comparison)
time nextflow run main.nf -profile test,docker --taxonomic_profiler kraken2 -resume

# New workflow
time nextflow run main.nf -profile test,docker --taxonomic_profiler kraken2 -resume
```

**Expected improvement**: 35-45% faster execution

#### Compare memory usage
```bash
# Check Nextflow execution report
cat results/pipeline_info/execution_report.html

# Look for:
# - Peak memory usage (should be 25-30% lower)
# - Process completion times (should be faster)
```

### 6. Compatibility Validation

#### Test with different databases
- Kraken2: SILVA, GTDB, RefSeq
- Sourmash: GTDB, GenBank

#### Test with different sample sizes
- Single sample
- Multiple samples (2-10)
- Large batch (>20 samples)

#### Test with different taxonomic levels
```bash
# Minimal levels
--taxonomy_plot_levels "Phylum,Genus"

# All levels
--taxonomy_plot_levels "Kingdom,Phylum,Class,Order,Family,Genus,Species"

# Custom selection
--taxonomy_plot_levels "Family,Genus,Species"
```

### 7. Error Handling Validation

#### Test with missing files
```bash
# Should fail gracefully with clear error message
python3 bin/taxonomy_report.py \
    --profiler kraken2 \
    --reports nonexistent.tsv \
    --reads-report Reads_report.csv \
    --db-name silva \
    --output-dir test_output/
```

#### Test with invalid profiler
```bash
# Should fail with validation error
python3 bin/taxonomy_report.py \
    --profiler invalid_profiler \
    --reports report.tsv \
    --reads-report Reads_report.csv \
    --db-name silva \
    --output-dir test_output/
```

#### Test with empty data
```bash
# Should handle gracefully
python3 bin/taxonomy_phyloseq.py \
    --profiler kraken2 \
    --input-files empty.biom \
    --db-name silva \
    --output-dir test_output/
```

## Success Criteria

### Functional Requirements
- ✅ All scripts execute without errors
- ✅ All modules complete successfully
- ✅ Output files are generated in correct locations
- ✅ Output formats match specifications

### Performance Requirements
- ✅ 35-45% faster execution time vs old workflow
- ✅ 25-30% lower memory usage vs old workflow
- ✅ No increase in disk I/O

### Quality Requirements
- ✅ TSV tables are valid and loadable in both Python and R
- ✅ HDF5 files are valid and readable
- ✅ R phyloseq objects (if generated) are valid
- ✅ Plots are generated and properly formatted
- ✅ No data loss or corruption

### Compatibility Requirements
- ✅ Works with both Kraken2 and Sourmash
- ✅ Compatible with existing pipeline parameters
- ✅ Backward compatible output formats (via optional R converter)
- ✅ Works across different execution profiles (docker, singularity, conda)

## Troubleshooting

### Common Issues

**Issue**: Python dependencies not found
```bash
# Solution: Install required packages
pip install pandas matplotlib numpy h5py seaborn
# Or use conda
conda install -c conda-forge pandas matplotlib numpy h5py seaborn
```

**Issue**: R phyloseq package not found
```bash
# Solution: Install phyloseq
R -e "install.packages('BiocManager'); BiocManager::install('phyloseq')"
```

**Issue**: Plots not generated
- Check that taxonomic levels exist in data
- Verify `--plot-levels` parameter is correct
- Check for sufficient classified taxa (need at least 1 taxon per level)

**Issue**: HDF5 file cannot be opened
- Verify h5py version compatibility
- Check file permissions
- Ensure sufficient disk space

## Reporting Issues

If validation fails, please report with:
1. Nextflow version
2. Execution profile used
3. Error messages (full stack trace)
4. Input data characteristics (sample count, database used)
5. Expected vs actual outputs

## References

- Optimization plan: `docs/taxonomy_optimization_plan.md`
- Migration guide: `docs/taxonomy_migration_guide.md`
- Pipeline manual: `docs/manual.md`
