# Phase 2: Output Structure Validation - COMPLETED

## Objective
Validate that assembly and binning module outputs follow the proper directory structure defined in `config/modules.config`.

## Audit Results

### ✅ Assembly Modules - COMPLIANT

#### MEGAHIT
- **Config Definition** (`config/modules.config:334-340`):
  ```groovy
  withName: 'MEGAHIT' {
      publishDir = [
          path: { meta.id == 'coassembly' ? 
              "${params.output}/03_assembly/coassembly" : 
              "${params.output}/03_assembly/per_sample/${meta.id}" },
          mode: params.publish_dir_mode,
          pattern: '*_contigs.fa'
      ]
  }
  ```

- **Module Output** (`modules/local/megahit/main.nf:33-34`):
  ```groovy
  tuple val(meta), path(reads), path("${meta.id}_contigs.fa"), emit: contigs_and_reads
  tuple val(meta), path("${meta.id}_contigs.fa"), emit: contigs_only
  ```

- **Status**: ✅ **CORRECT** - Outputs `*_contigs.fa` matching pattern

#### BBMAP
- **Config Definition** (`config/modules.config:343-356`):
  ```groovy
  withName: 'BBMAP' {
      publishDir = [
          [
              path: { meta.id == 'coassembly' ? 
                  "${params.output}/03_assembly/coassembly" : 
                  "${params.output}/03_assembly/per_sample/${meta.id}" },
              mode: params.publish_dir_mode,
              pattern: '*_filtered_contigs.fa'
          ],
          [
              path: { meta.id == 'coassembly' ? 
                  "${params.output}/03_assembly/coassembly" : 
                  "${params.output}/03_assembly/per_sample/${meta.id}" },
              mode: params.publish_dir_mode,
              pattern: '*_contig.stats'
          ]
      ]
  }
  ```

- **Module Output** (`modules/local/bbmap/main.nf:32-36`):
  ```groovy
  tuple val(meta), path(reads), path("${meta.id}_filtered_contigs.fa"), emit: contigs_with_reads
  tuple val(meta), path("${meta.id}_filtered_contigs.fa"), emit: contigs_only
  path "${meta.id}_contig.stats", emit: stats
  path "${meta.id}_filter_report.txt", emit: filter_report
  ```

- **Status**: ✅ **CORRECT** - Outputs match patterns

### ✅ Binning Modules - COMPLIANT

#### METABAT2
- **Config Definition** (`config/modules.config:363-371`):
  ```groovy
  withName: 'METABAT2' {
      publishDir = [
          path: { meta.id == 'coassembly' ? 
              "${params.output}/04_binning/coassembly/raw_bins/metabat2" : 
              "${params.output}/04_binning/per_sample/${meta.id}/raw_bins/metabat2" },
          mode: params.publish_dir_mode,
          pattern: '*metabat_bins'
      ]
  }
  ```

- **Module Output** (`modules/local/metabat2/main.nf:30`):
  ```groovy
  tuple val(meta), path("${meta.id}_metabat_bins"), emit: bins
  ```

- **Status**: ✅ **CORRECT** - Outputs `*metabat_bins` matching pattern

#### SEMIBIN
- **Config Definition** (`config/modules.config:373-381`):
  ```groovy
  withName: 'SEMIBIN' {
      publishDir = [
          path: { meta.id == 'coassembly' ? 
              "${params.output}/04_binning/coassembly/raw_bins/semibin" : 
              "${params.output}/04_binning/per_sample/${meta.id}/raw_bins/semibin" },
          mode: params.publish_dir_mode,
          pattern: '*semibin_output_bins'
      ]
  }
  ```

- **Module Output** (`modules/local/semibin/main.nf:32`):
  ```groovy
  tuple val(meta), path("${meta.id}_semibin_output_bins"), emit: bins
  ```

- **Status**: ✅ **CORRECT** - Outputs `*semibin_output_bins` matching pattern

#### COMEBIN
- **Config Definition** (`config/modules.config:383-391`):
  ```groovy
  withName: 'COMEBIN' {
      publishDir = [
          path: { meta.id == 'coassembly' ? 
              "${params.output}/04_binning/coassembly/raw_bins/comebin" : 
              "${params.output}/04_binning/per_sample/${meta.id}/raw_bins/comebin" },
          mode: params.publish_dir_mode,
          pattern: '*_comebin_bins/comebin_res/comebin_res_bins'
      ]
  }
  ```

- **Module Output** (`modules/local/comebin/main.nf:31`):
  ```groovy
  tuple val(meta), path("${meta.id}_comebin_bins/comebin_res/comebin_res_bins"), emit: bins
  ```

- **Status**: ✅ **CORRECT** - Outputs match pattern

#### METAWRAP
- **Config Definition** (`config/modules.config:393-401`):
  ```groovy
  withName: 'METAWRAP' {
      publishDir = [
          path: { meta.id == 'coassembly' ? 
              "${params.output}/04_binning/coassembly/refined_bins" : 
              "${params.output}/04_binning/per_sample/${meta.id}/refined_bins" },
          mode: params.publish_dir_mode,
          pattern: '*metawrap*bins'
      ]
  }
  ```

- **Module Output** (`modules/local/metawrap/main.nf:31`):
  ```groovy
  tuple val(meta), path("${meta.id}_metawrap_${params.metawrap_completeness}_${params.metawrap_contamination}_bins"), emit: bins
  ```

- **Status**: ✅ **CORRECT** - Outputs `*metawrap*bins` matching pattern

#### CHECKM2
- **Config Definition** (`config/modules.config:403-411`):
  ```groovy
  withName: 'CHECKM2' {
      publishDir = [
          path: { meta.id == 'coassembly' ? 
              "${params.output}/04_binning/coassembly/quality/checkm2" : 
              "${params.output}/04_binning/per_sample/${meta.id}/quality/checkm2" },
          mode: params.publish_dir_mode,
          pattern: '*quality_report.tsv'
      ]
  }
  ```

- **Module Output** (`modules/local/checkm2/main.nf:31-32`):
  ```groovy
  tuple val(meta), path("*quality_report.tsv"), emit: all_reports
  tuple val(meta), path("*_metawrap_quality_report.tsv"), emit: metawrap_report
  ```

- **Status**: ✅ **CORRECT** - Outputs match pattern

#### GTDB_TK
- **Config Definition** (`config/modules.config:413-421`):
  ```groovy
  withName: 'GTDB_TK' {
      publishDir = [
          path: { meta.id == 'coassembly' ? 
              "${params.output}/04_binning/coassembly/taxonomy/gtdbtk" : 
              "${params.output}/04_binning/per_sample/${meta.id}/taxonomy/gtdbtk" },
          mode: params.publish_dir_mode,
          pattern: '*_gtdbtk_*'
      ]
  }
  ```

- **Module Output** (`modules/local/gtdb-tk/main.nf:33-34`):
  ```groovy
  tuple val(meta), path("*_gtdbtk_*"), emit: gtdb_tk
  tuple val(meta), path("*_gtdbtk_*"), emit: report
  ```

- **Status**: ✅ **CORRECT** - Outputs match pattern

### ✅ Report Modules - COMPLIANT

#### BIN_QUALITY_REPORT
- **Config**: `04_binning/quality/summary` for `*.png` and `*.csv`
- **Status**: ✅ **CORRECT** - Config properly defined

#### BIN_TAX_REPORT
- **Config**: `04_binning/taxonomy/summary` for `*.png` and `*.csv`
- **Status**: ✅ **CORRECT** - Config properly defined

#### BIN_SUMMARY
- **Config**: `04_binning/coassembly/summary` for `Bin_summary.csv`
- **Status**: ✅ **CORRECT** - Config properly defined

## Expected Output Directory Structure

```
{output}/
├── 01_quality_control/
│   ├── fastp/{sample_id}/
│   └── summary/
├── 02_taxonomy/
│   ├── kraken2/
│   ├── sourmash/
│   ├── figures/
│   ├── tables/
│   └── phyloseq/
├── 03_assembly/
│   ├── per_sample/{sample_id}/
│   │   ├── {sample_id}_contigs.fa
│   │   ├── {sample_id}_filtered_contigs.fa
│   │   └── {sample_id}_contig.stats
│   └── coassembly/
│       ├── coassembly_contigs.fa
│       ├── coassembly_filtered_contigs.fa
│       └── coassembly_contig.stats
├── 04_binning/
│   ├── per_sample/{sample_id}/
│   │   ├── raw_bins/
│   │   │   ├── metabat2/
│   │   │   ├── semibin/
│   │   │   └── comebin/
│   │   ├── refined_bins/
│   │   ├── quality/checkm2/
│   │   └── taxonomy/gtdbtk/
│   ├── coassembly/
│   │   ├── raw_bins/
│   │   │   ├── metabat2/
│   │   │   ├── semibin/
│   │   │   └── comebin/
│   │   ├── refined_bins/
│   │   ├── quality/checkm2/
│   │   ├── taxonomy/gtdbtk/
│   │   ├── coverage/
│   │   └── summary/
│   ├── quality/summary/
│   └── taxonomy/summary/
├── 05_arg_prediction/
├── 06_contig_taxonomy/
├── 07_functional_annotation/
└── pipeline_info/
```

## Summary

✅ **All assembly and binning modules are correctly configured**

- All module outputs match their corresponding publishDir patterns in `config/modules.config`
- Both per-sample and co-assembly modes use proper conditional paths based on `meta.id`
- Output file naming conventions are consistent
- No fixes required for output structure

## Observations

1. **Unified Module Design**: The refactored modules properly handle both per-sample and co-assembly modes using `meta.id == 'coassembly'` checks
2. **Pattern Matching**: All publishDir patterns correctly match the actual output file names
3. **Directory Hierarchy**: The config follows a logical structure separating raw bins, refined bins, quality reports, and taxonomy
4. **Consistency**: Naming conventions are consistent across all modules

## Issues Found and Fixed

After initial audit, user correctly identified additional unwanted output folders being created:
- `bowtie2/` - from BOWTIE2_SAMTOOLS processes
- `calculate/` - from CALCULATE_DEPTH process

### Root Cause
Default publishDir in `conf/base.config` applies to ALL processes unless explicitly overridden in `config/modules.config`.

### Processes Fixed
Added publishDir configurations for 7 processes that were missing explicit configs:

1. **COUNT_READS** - Disabled (intermediate, counts in READS_REPORT)
2. **BOWTIE2_SAMTOOLS** - Disabled (intermediate BAM files)
3. **BOWTIE2_SAMTOOLS_DEPTH** - Disabled (intermediate BAM files)
4. **CONTIG_FILTER_SUMMARY** - Disabled (summary in other reports)
5. **CALCULATE_DEPTH** - Disabled (intermediate depth files)
6. **NT_BLASTN** - Disabled (intermediate BLAST results)
7. **BLOBTOOLS** - Disabled (intermediate blob tables)

See `PHASE2_FIXES_APPLIED.md` for detailed changes.

## Status
✅ **PHASE 2 COMPLETED** - All output structure issues resolved. All modules now have proper publishDir configurations.

## Next Steps
Proceed to Phase 3: Custom container migration to stable containers with Python scripts.
