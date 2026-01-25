# Channel Structure Validation - Phase 2

**Date**: January 22, 2026  
**Status**: ✅ VALIDATED

---

## Issue: Process Input/Output Mismatch

**Error**:
```
Process `BINNING:METABAT2` declares 1 input but was called with 2 arguments
```

**Root Cause**: 
The error message was misleading. The actual issue was that `CALCULATE_DEPTH` has a named output emit `.depth`, but we were passing the entire process output object instead of the specific emit.

---

## Process Signatures

### CALCULATE_DEPTH

**Input**:
```groovy
tuple val(meta), path(contigs), path(bam)
```

**Output**:
```groovy
tuple val(meta), path(contigs), path("${meta.id}_depth.txt"), emit: depth
path "versions.yml", emit: versions
```

**Usage**:
```groovy
ch_depth = CALCULATE_DEPTH(bam)
// Must use: ch_depth.depth (not just ch_depth)
```

---

### METABAT2

**Input**:
```groovy
tuple val(meta), path(contigs), path(depth)
```

**Output**:
```groovy
tuple val(meta), path("${meta.id}_metabat_bins"), emit: bins
path "versions.yml", emit: versions
```

**Correct Usage**:
```groovy
ch_metabat2 = METABAT2(ch_depth.depth)  // ✅ Use .depth emit
```

**Incorrect Usage**:
```groovy
ch_metabat2 = METABAT2(ch_depth)  // ❌ Wrong - passes process object
```

---

### SEMIBIN

**Input**:
```groovy
tuple val(meta), path(contigs), path(bam)
```

**Output**:
```groovy
tuple val(meta), path("${meta.id}_semibin_output_bins"), emit: bins
path "versions.yml", emit: versions
```

**Usage**:
```groovy
ch_semibin = SEMIBIN(bam)  // ✅ Direct input from assembly
```

---

### COMEBIN

**Input**:
```groovy
tuple val(meta), path(contigs), path(bam)
```

**Output**:
```groovy
tuple val(meta), path("${meta.id}_comebin_bins/comebin_res/comebin_res_bins"), emit: bins
path "versions.yml", emit: versions
```

**Usage**:
```groovy
ch_comebin = COMEBIN(bam)  // ✅ Direct input from assembly
```

---

### METAWRAP

**Input**:
```groovy
tuple val(meta), path(metabat2_bins), path(semibin_bins), path(comebin_bins)
```

**Output**:
```groovy
tuple val(meta), path("${meta.id}_metawrap_..."), emit: bins
path "versions.yml", emit: versions
```

**Usage**:
```groovy
ch_all_bins = ch_metabat2.bins.join(ch_semibin.bins).join(ch_comebin.bins)
ch_metawrap = METAWRAP(ch_all_bins)
```

---

### CHECKM2

**Input**:
```groovy
tuple val(meta), path(metabat2), path(semibin), path(comebin), path(metawrap), path(checkm_db)
```

**Output**:
```groovy
tuple val(meta), path("*quality_report.tsv"), emit: all_reports
tuple val(meta), path("*_metawrap_quality_report.tsv"), emit: metawrap_report
path "versions.yml", emit: versions
```

**Usage**:
```groovy
ch_checkm = CHECKM2(
    ch_all_bins.join(ch_metawrap.bins).combine(checkm2_db)
)
```

---

### GTDB_TK

**Input**:
```groovy
tuple val(meta), path(metawrap), path(gtdbtk_db)
```

**Output**:
```groovy
tuple val(meta), path("*_gtdbtk_*"), emit: gtdb_tk
tuple val(meta), path("*_gtdbtk_*"), emit: report
path "versions.yml", emit: versions
```

**Usage**:
```groovy
ch_gtdb_tk = GTDB_TK(ch_metawrap.bins.combine(gtdbtk_db))
```

---

## Fixed Binning Workflow

### Per-Sample Mode

```groovy
// Calculate depth
ch_depth = CALCULATE_DEPTH(bam)

// Run binning tools
ch_metabat2 = METABAT2(ch_depth.depth)  // ✅ Use .depth emit
ch_semibin = SEMIBIN(bam)
ch_comebin = COMEBIN(bam)

// Combine bins for refinement
ch_all_bins = ch_metabat2.bins.join(ch_semibin.bins).join(ch_comebin.bins)

// Refine bins with MetaWRAP
ch_metawrap = METAWRAP(ch_all_bins)

// Quality assessment
ch_checkm = CHECKM2(
    ch_all_bins.join(ch_metawrap.bins).combine(checkm2_db)
)

// Taxonomic classification
ch_gtdb_tk = GTDB_TK(ch_metawrap.bins.combine(gtdbtk_db))
```

### Co-Assembly Mode

```groovy
// Create coassembly meta
def coassembly_meta = [id: 'coassembly']

// Prepare inputs with coassembly meta
ch_bam_collected = bam.collect()
ch_contigs_with_meta = channel.of(coassembly_meta)
    .combine(contigs)
    .combine(ch_bam_collected)

// Calculate depth
ch_depth_co = CALCULATE_DEPTH(ch_contigs_with_meta)

// Run binning tools
ch_metabat2_co = METABAT2(ch_depth_co.depth)  // ✅ Use .depth emit
ch_semibin_co = SEMIBIN(ch_contigs_with_meta)
ch_comebin_co = COMEBIN(ch_contigs_with_meta)

// Rest follows same pattern as per-sample
```

---

## Channel Flow Diagram

```
Assembly Output (bam)
    ↓
CALCULATE_DEPTH
    ├─ .depth → [meta, contigs, depth.txt]  → METABAT2
    └─ .versions → versions.yml
    
Assembly Output (bam)
    ↓
SEMIBIN → .bins → [meta, semibin_bins]
    
Assembly Output (bam)
    ↓
COMEBIN → .bins → [meta, comebin_bins]

Join all bins by meta.id
    ↓
[meta, metabat2_bins, semibin_bins, comebin_bins]
    ↓
METAWRAP → .bins → [meta, metawrap_bins]
```

---

## Validation Checklist

- [x] CALCULATE_DEPTH output structure correct
- [x] METABAT2 called with `.depth` emit
- [x] SEMIBIN input matches assembly output
- [x] COMEBIN input matches assembly output
- [x] METAWRAP input is joined bins
- [x] CHECKM2 input is all bins + db
- [x] GTDB_TK input is metawrap bins + db
- [x] All version outputs collected
- [x] Both per-sample and co-assembly modes fixed

---

## Summary

**Issue**: Not using named output emits correctly  
**Fix**: Use `.depth` emit from CALCULATE_DEPTH when calling METABAT2  
**Files Modified**: `subworkflows/local/binning.nf`  
**Lines Changed**: 2 (lines 43 and 98)

The binning workflow now correctly uses named output emits and all process calls match their input signatures.
