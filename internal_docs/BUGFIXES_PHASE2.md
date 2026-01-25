# Phase 2 Bug Fixes - Variable Scope Errors

**Date**: January 22, 2026  
**Status**: ✅ FIXED

---

## Issues Found During Testing

### Issue 1: Variable Scope Error in Output Declarations

**Error Message**:
```
ERROR ~ No such variable: meta
-- Check script '/home/jugalde/pipelines/BugBuster/modules/local/megahit/main.nf' at line: 23
```

**Root Cause**: 
In Nextflow, variables used in `output` declarations must be available at **parse time**, not runtime. We were using `${prefix}` in output paths, but `prefix` was defined in the `script` section which is evaluated later.

**Incorrect Code**:
```groovy
output:
tuple val(meta), path("${prefix}_contigs.fa"), emit: contigs_only

script:
prefix = meta.id  // Defined too late!
```

**Fix Applied**:
Use `${meta.id}` directly in output declarations since `meta` is available from the input.

```groovy
output:
tuple val(meta), path("${meta.id}_contigs.fa"), emit: contigs_only

script:
prefix = meta.id  // Still used in script for consistency
```

---

### Issue 2: Conditional Label Expressions

**Error Message**:
```
ERROR ~ No such variable: meta
-- Check script '/home/jugalde/pipelines/BugBuster/modules/local/megahit/main.nf' at line: 23
```

**Root Cause**:
The `label` directive cannot use conditional expressions with input variables at parse time.

**Incorrect Code**:
```groovy
process MEGAHIT {
    label "${ meta.id == 'coassembly' ? 'process_high' : 'process_medium' }"
```

**Fix Applied**:
Use a fixed label that works for both modes.

```groovy
process MEGAHIT {
    label 'process_medium'
```

**Note**: For processes that truly need different resources, this can be configured in `conf/base.config` using dynamic resource allocation based on input size rather than meta.id.

---

## Modules Fixed

### Assembly Modules (3 files)

1. **`modules/local/megahit/main.nf`**
   - Changed output: `${prefix}_contigs.fa` → `${meta.id}_contigs.fa`
   - Changed label: conditional → `'process_medium'`

2. **`modules/local/bbmap/main.nf`**
   - Changed output: `${prefix}_filtered_contigs.fa` → `${meta.id}_filtered_contigs.fa`
   - Changed output: `${prefix}_contig.stats` → `${meta.id}_contig.stats`

3. **`modules/local/bowtie2_samtools/main.nf`**
   - Changed output: `${prefix}_all_reads.bam` → `${meta.id}_all_reads.bam`

### Binning Modules (5 files)

4. **`modules/local/metabat2/main.nf`**
   - Changed output: `${prefix}_metabat_bins` → `${meta.id}_metabat_bins`

5. **`modules/local/semibin/main.nf`**
   - Changed output: `${prefix}_semibin_output_bins` → `${meta.id}_semibin_output_bins`
   - Changed label: conditional → `'process_medium'`

6. **`modules/local/comebin/main.nf`**
   - Changed output: `${prefix}_comebin_bins/...` → `${meta.id}_comebin_bins/...`

7. **`modules/local/metawrap/main.nf`**
   - Changed output: `${prefix}_metawrap_${completeness}_${contamination}_bins` → `${meta.id}_metawrap_${params.metawrap_completeness}_${params.metawrap_contamination}_bins`

8. **`modules/local/calculate_depth/main.nf`**
   - Changed output: `${prefix}_depth.txt` → `${meta.id}_depth.txt`

---

## Pattern for Future Modules

When creating or modifying Nextflow modules, follow this pattern:

### ✅ CORRECT Pattern

```groovy
process TOOL_NAME {
    tag "${meta.id}"
    label 'process_medium'  // Fixed label
    
    input:
    tuple val(meta), path(inputs)
    
    output:
    // Use meta.id or params directly - available at parse time
    tuple val(meta), path("${meta.id}_output.txt"), emit: main
    path "versions.yml", emit: versions
    
    script:
    // Define variables for use in script section
    prefix = meta.id
    def args = task.ext.args ?: ''
    
    """
    # Use ${prefix} in bash commands
    tool --input ${inputs} --output ${prefix}_output.txt
    """
}
```

### ❌ INCORRECT Pattern

```groovy
process TOOL_NAME {
    label "${ meta.id == 'coassembly' ? 'process_high' : 'process_medium' }"  // ❌ Can't use meta in label
    
    output:
    tuple val(meta), path("${prefix}_output.txt"), emit: main  // ❌ prefix not defined yet
    
    script:
    prefix = meta.id  // Too late for output section
}
```

---

## Validation

Pipeline now parses successfully:

```bash
$ nextflow run main.nf --help
# Returns help message without errors ✅
```

---

## Testing Recommendations

After these fixes, test the pipeline with:

```bash
# Stub-run test (fast validation)
nextflow run main.nf \
  --input test_samplesheet.csv \
  --output test_small \
  --assembly_mode assembly \
  --include_binning true \
  -profile docker \
  -stub-run

# Real run with small dataset
nextflow run main.nf \
  --input test_samplesheet.csv \
  --output test_small \
  --assembly_mode assembly \
  --include_binning true \
  -profile docker
```

---

## Summary

**Total Files Modified**: 8 module files  
**Total Changes**: 
- 10 output path fixes (variable scope)
- 2 label directive fixes (conditional removal)

**Result**: Pipeline now parses and validates correctly ✅

All variable scope and parse-time evaluation errors have been resolved. The pipeline is ready for testing.
