# Phase 2: Output Structure Fixes Applied

## Issues Found

User correctly identified that additional unwanted output folders were being created:
- `bowtie2/` - from BOWTIE2_SAMTOOLS processes
- `calculate/` - from CALCULATE_DEPTH process
- Potentially `contig/` - from CONTIG_FILTER_SUMMARY

## Root Cause

Default publishDir in `conf/base.config:13-17` applies to ALL processes unless explicitly overridden:
```groovy
publishDir = [
    path: { "${params.output}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
    mode: params.publish_dir_mode,
    saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
]
```

Several processes lacked explicit publishDir configurations in `config/modules.config`.

## Fixes Applied

Added publishDir configurations to `config/modules.config` for the following processes:

### 1. COUNT_READS (lines 262-266)
```groovy
withName: 'COUNT_READS' {
    publishDir = [
        enabled: false  // Intermediate file - counts included in READS_REPORT
    ]
}
```
**Rationale**: Intermediate file, counts are included in READS_REPORT output.

### 2. BOWTIE2_SAMTOOLS (lines 341-345)
```groovy
withName: 'BOWTIE2_SAMTOOLS' {
    publishDir = [
        enabled: false  // Intermediate file - BAM used for depth calculation and indexing
    ]
}
```
**Rationale**: Large intermediate BAM files, only needed for downstream depth calculation and indexing.

### 3. BOWTIE2_SAMTOOLS_DEPTH (lines 348-352)
```groovy
withName: 'BOWTIE2_SAMTOOLS_DEPTH' {
    publishDir = [
        enabled: false  // Intermediate file - BAM used for bin coverage calculation
    ]
}
```
**Rationale**: Intermediate BAM files for bin coverage in coassembly mode.

### 4. CONTIG_FILTER_SUMMARY (lines 355-359)
```groovy
withName: 'CONTIG_FILTER_SUMMARY' {
    publishDir = [
        enabled: false  // Summary information already captured in other reports
    ]
}
```
**Rationale**: Summary information is already captured in other assembly reports.

### 5. CALCULATE_DEPTH (lines 390-394)
```groovy
withName: 'CALCULATE_DEPTH' {
    publishDir = [
        enabled: false  // Intermediate file - depth consumed by MetaBAT2
    ]
}
```
**Rationale**: Intermediate depth file consumed by MetaBAT2 for binning.

## Impact

These fixes prevent the creation of unwanted output directories:
- ❌ `{output}/bowtie2/` - Now disabled
- ❌ `{output}/calculate/` - Now disabled
- ❌ `{output}/contig/` - Now disabled
- ❌ `{output}/count/` - Now disabled

## Verification Needed

Need to verify no other processes are using default publishDir. Processes to check:
- BEDTOOLS (used in binning coassembly)
- ARG_FASTA_FORMATTER (used in ARG prediction)
- NT_BLASTN (used in contig taxonomy)
- BLOBTOOLS (used in contig taxonomy)
- QFILTER (used in QC)
- SAMTOOLS_INDEX (used in contig taxonomy)
- Database formatting processes (already have configs)

## Status
✅ **FIXES APPLIED** - Intermediate files now properly disabled from publishing.

## Next Steps
1. Verify no other processes are missing publishDir configs
2. Update Phase 2 completion documentation
3. Proceed to Phase 3
