# Phase 2: Output Structure Issues - CRITICAL FINDINGS

## Problem Identified

The user correctly identified that additional output folders (`bowtie2` and `calculate`) are being created that were not documented in the initial Phase 2 audit.

## Root Cause

**Default publishDir in `conf/base.config:13-17`:**
```groovy
process {
    // Default publish mode
    publishDir = [
        path: { "${params.output}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

This default publishDir applies to **ALL processes** unless explicitly overridden in `config/modules.config`.

## Missing publishDir Configurations

The following processes are **NOT** configured in `config/modules.config` and are falling back to the default:

### 1. CALCULATE_DEPTH
- **Process**: `CALCULATE_DEPTH`
- **Default Output Path**: `${params.output}/calculate/`
- **Files Created**: `{sample_id}_depth.txt`
- **Used By**: Binning subworkflow (per-sample mode)
- **Should Output To**: Should be disabled (intermediate file) OR `04_binning/per_sample/{sample_id}/depth/`

### 2. BOWTIE2_SAMTOOLS
- **Process**: `BOWTIE2_SAMTOOLS`
- **Default Output Path**: `${params.output}/bowtie2/`
- **Files Created**: `{sample_id}_all_reads.bam`
- **Used By**: Assembly subworkflow (for binning and contig analysis)
- **Should Output To**: Should be disabled (intermediate file) OR `03_assembly/per_sample/{sample_id}/alignment/`

### 3. BOWTIE2_SAMTOOLS_DEPTH
- **Process**: `BOWTIE2_SAMTOOLS_DEPTH`
- **Default Output Path**: `${params.output}/bowtie2/`
- **Files Created**: `{sample_id}_{bin_name}_all_reads.bam` (multiple BAMs per bin)
- **Used By**: Binning subworkflow (co-assembly mode for bin coverage)
- **Should Output To**: Should be disabled (intermediate file)

### 4. CONTIG_FILTER_SUMMARY
- **Process**: `CONTIG_FILTER_SUMMARY`
- **Default Output Path**: `${params.output}/contig/`
- **Files Created**: Summary reports
- **Used By**: Assembly subworkflow
- **Should Output To**: `03_assembly/summary/` OR disabled

## Additional Processes Without Explicit Config

Checking all local modules, the following may also be using default publishDir:

- **COUNT_READS** - Used in QC subworkflow (should be disabled - intermediate)
- **CONTIG_FILTER_SUMMARY** - Used in assembly subworkflow
- **BEDTOOLS** - Used in binning subworkflow (co-assembly)
- **ARG_FASTA_FORMATTER** - Used in ARG prediction
- **COMPUTE_KMER_FREQ** - If used anywhere
- **BASALT** - If used anywhere
- **AUTOMETA** - If used anywhere
- **VAMB** - If used anywhere

## Impact

These processes are creating unwanted output directories:
```
{output}/
├── bowtie2/          ← UNWANTED (from BOWTIE2_SAMTOOLS)
├── calculate/        ← UNWANTED (from CALCULATE_DEPTH)
├── contig/           ← UNWANTED (from CONTIG_FILTER_SUMMARY)
└── ... (potentially more)
```

## Recommended Fixes

### Option 1: Disable Publishing (Recommended for Intermediate Files)
Add to `config/modules.config`:

```groovy
withName: 'CALCULATE_DEPTH' {
    publishDir = [
        enabled: false  // Intermediate file - depth used by MetaBAT2
    ]
}

withName: 'BOWTIE2_SAMTOOLS' {
    publishDir = [
        enabled: false  // Intermediate file - BAM used for depth/indexing
    ]
}

withName: 'BOWTIE2_SAMTOOLS_DEPTH' {
    publishDir = [
        enabled: false  // Intermediate file - BAM used for bin coverage
    ]
}

withName: 'CONTIG_FILTER_SUMMARY' {
    publishDir = [
        enabled: false  // Summary already captured elsewhere
    ]
}

withName: 'COUNT_READS' {
    publishDir = [
        enabled: false  // Intermediate file - counts in reports
    ]
}
```

### Option 2: Proper Directory Structure (If Files Should Be Published)
Add to `config/modules.config`:

```groovy
withName: 'CALCULATE_DEPTH' {
    publishDir = [
        path: { meta.id == 'coassembly' ? 
            "${params.output}/04_binning/coassembly/depth" : 
            "${params.output}/04_binning/per_sample/${meta.id}/depth" },
        mode: params.publish_dir_mode,
        pattern: '*_depth.txt'
    ]
}

withName: 'BOWTIE2_SAMTOOLS' {
    publishDir = [
        path: { meta.id == 'coassembly' ? 
            "${params.output}/03_assembly/coassembly/alignment" : 
            "${params.output}/03_assembly/per_sample/${meta.id}/alignment" },
        mode: params.publish_dir_mode,
        pattern: '*.bam'
    ]
}

withName: 'CONTIG_FILTER_SUMMARY' {
    publishDir = [
        path: { "${params.output}/03_assembly/summary" },
        mode: params.publish_dir_mode,
        pattern: '*.{txt,csv}'
    ]
}
```

## Recommendation

**Disable publishing for all intermediate files** (Option 1) because:
1. BAM files are large and intermediate - only needed for downstream processes
2. Depth files are intermediate - consumed by MetaBAT2
3. These files clutter the output directory
4. Users can always access them in the work directory if needed

Only publish final results and reports.

## Status

❌ **PHASE 2 INCOMPLETE** - Additional fixes required before proceeding to Phase 3.

## Next Steps

1. Add publishDir configurations for all missing processes
2. Verify no other processes are using default publishDir
3. Re-validate complete output structure
4. Update PHASE2_COMPLETION.md with corrected information
