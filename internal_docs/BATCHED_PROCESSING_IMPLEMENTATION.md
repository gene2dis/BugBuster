# Batched CheckM2 and GTDB-Tk Implementation

## Overview

This document describes the implementation of batched processing for CheckM2 and GTDB-Tk modules to improve pipeline performance by reducing database loading overhead.

## Performance Benefits

### Database Loading Optimization
- **CheckM2**: Database (~3-4GB) loaded once instead of N times
- **GTDB-Tk**: Database (~80GB uncompressed) loaded once instead of N times
- **Estimated savings**:
  - 10 samples: 20-50 minutes
  - 50 samples: 2-4 hours
  - 100 samples: 4-8 hours

### Additional Benefits
- Reduced tool initialization overhead
- Better resource utilization (single large job vs multiple small jobs)
- Removed `maxForks 2` limitation on GTDB-Tk
- Improved parallelization efficiency

## Implementation Details

### New Modules Created

#### 1. CHECKM2_BATCH (`modules/local/checkm2_batch/main.nf`)
- **Input**: `tuple val(meta_list), path(all_bins), path(checkm_db)`
- **Process**:
  - Collects bins from all samples with prefixed names (`${sample_id}__${bin_name}.fa`)
  - Runs CheckM2 once per binner across all samples
  - Splits combined output back to per-sample reports
- **Output**: Per-sample quality reports matching original format

#### 2. GTDB_TK_BATCH (`modules/local/gtdb_tk_batch/main.nf`)
- **Input**: `tuple val(meta_list), path(all_bins), path(gtdbtk_db)`
- **Process**:
  - Collects metawrap bins from all samples with prefixed names
  - Runs GTDB-Tk once across all samples
  - Splits combined output (bac120 and ar53) back to per-sample reports
- **Output**: Per-sample taxonomy reports matching original format

### Modified Files

#### 1. Binning Subworkflow (`subworkflows/local/binning.nf`)

**Changes**:
- Replaced `CHECKM2` with `CHECKM2_BATCH`
- Replaced `GTDB_TK` with `GTDB_TK_BATCH`
- Added channel transformations to:
  - Collect all bins from all samples
  - Process in batch
  - Split results back to per-sample format using `flatMap`

**Channel Flow**:
```groovy
// CheckM2 batching
ch_all_sample_bins = ch_metabat2.bins
    .join(ch_semibin.bins)
    .join(ch_comebin.bins)
    .join(ch_metawrap.bins)
    .collect()  // Collect all samples
    .map { items -> [meta_list, all_bins] }  // Organize for batch processing

ch_checkm_batch = CHECKM2_BATCH(ch_all_sample_bins.combine(checkm2_db))

// Transform back to per-sample
ch_checkm_all_reports = ch_checkm_batch.all_reports
    .flatMap { meta_list, reports ->
        // Split reports by sample prefix
        sample_reports.collect { sample_id, report_list ->
            [meta, report_list]
        }
    }
```

## Output Structure Preservation

### Key Design Principle
Results are processed in batch but split back to maintain the original per-sample/coassembly output structure.

### Bin Naming Convention
- **During processing**: `${sample_id}__${original_bin_name}.fa`
- **In output**: `${original_bin_name}.fa` (prefix removed)

### Report Splitting
Python scripts embedded in the modules parse combined reports and split by sample prefix:
- Extract sample ID from bin name
- Group results by sample
- Write individual per-sample reports with original bin names

## Compatibility

### Downstream Modules
No changes required to downstream modules:
- `BIN_QUALITY_REPORT`
- `BIN_TAX_REPORT`
- `BIN_SUMMARY`

All receive the same channel structure and file formats as before.

### Assembly Modes
Works seamlessly with both:
- **Assembly mode**: Per-sample processing
- **Co-assembly mode**: Co-assembly processing

## Testing Recommendations

### Test Scenarios
1. **Single sample**: Verify no regression
2. **Multiple samples (5-10)**: Verify batching works correctly
3. **Co-assembly mode**: Verify compatibility
4. **Empty bins**: Verify empty report handling
5. **Mixed bacteria/archaea**: Verify both bac120 and ar53 outputs

### Validation Checks
- Output format matches original exactly
- All sample reports generated
- Downstream reports work without modification
- Performance improvement measured
- Resource usage monitored

### Performance Testing
```bash
# Run with timing
time nextflow run main.nf -profile test --include_binning

# Compare with previous implementation
# Expected: Significant time reduction for multiple samples
```

## Troubleshooting

### Common Issues

1. **Missing reports for some samples**
   - Check bin naming convention
   - Verify sample ID extraction logic in Python scripts

2. **Empty reports**
   - Normal behavior when no bins found
   - Matches original module behavior

3. **Memory issues**
   - GTDB-Tk requires high memory (labeled `process_high_memory`)
   - CheckM2 labeled `process_medium`

## Future Enhancements

### Optional Configuration
Could add parameter to toggle batching:
```groovy
params.batch_checkm2_gtdbtk = true  // Set to false for per-sample processing
```

This would allow users to revert to per-sample processing if needed for debugging or specific use cases.

## Summary

The batched implementation provides significant performance improvements while maintaining complete compatibility with existing pipeline structure and downstream modules. The approach is transparent to end users and requires no changes to pipeline invocation or configuration.
