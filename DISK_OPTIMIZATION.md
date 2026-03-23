# BugBuster Pipeline - Disk Space Optimization

This document describes the disk space optimization features implemented in BugBuster to prevent running out of disk space during pipeline execution.

## Overview

The pipeline implements a **progressive cleanup strategy** that frees disk space **during execution**, not after. This is achieved through three complementary mechanisms:

1. **publishDir with mode 'move'** - Move (not copy) large outputs to final location, immediately freeing work directory space
2. **Internal process cleanup** - Remove temporary files within processes as soon as they're no longer needed
3. **Nextflow cleanup configuration** - Enable automatic work directory cleanup with resume support

## Implementation Details

### Phase 1: publishDir with mode 'move' for Large Intermediates

The following processes use `publishDir` with `mode: 'move'` to immediately free disk space by moving (not copying) outputs:

#### BOWTIE2 (Clean Reads)
- **Location**: `modules/local/bowtie2/main.nf`
- **Moves to**: `output/clean_reads/{sample_id}/`
- **When**: After phiX removal (final QC step)
- **Mode**: `move` - files are moved from work directory, not copied
- **Disk savings**: ~8-10 GB per sample
- **Files moved**: `*_R1_map_phiX.fastq.gz`, `*_R2_map_phiX.fastq.gz`, `*_bowtie_report.tsv`

#### BBMAP (Filtered Contigs)
- **Location**: `modules/local/bbmap/main.nf`
- **Moves to**: `output/assembly/{sample_id}/`
- **When**: After contig length filtering
- **Mode**: `move` - files are moved from work directory, not copied
- **Disk savings**: ~2-5 GB per sample
- **Files moved**: `*_filtered_contigs.fa`, `*_contig.stats`, `*_filter_report.txt`

#### METAWRAP (Refined Bins)
- **Location**: `modules/local/metawrap/main.nf`
- **Moves to**: `output/bins/{sample_id}/`
- **When**: After bin refinement
- **Mode**: `move` - files are moved from work directory, not copied
- **Disk savings**: ~1-3 GB per sample
- **Files moved**: `*_metawrap_*_bins/` directory with all refined bins

### Phase 2: Internal Process Cleanup

The following processes clean up temporary files during execution:

#### MEGAHIT (Assembly)
- **Location**: `modules/local/megahit/main.nf`
- **Removes**:
  - `intermediate_contigs/` - Intermediate assembly files
  - `kmer_k*/` - K-mer data directories
  - `*.tmp` - Temporary files
  - `checkpoints.txt` - Assembly checkpoints
- **Keeps**: Final contigs and log files
- **Disk savings**: ~1-3 GB per sample
- **Timing**: Immediately after moving final contigs

#### BOWTIE2_SAMTOOLS (Alignment)
- **Location**: `modules/local/bowtie2_samtools/main.nf`
- **Removes**:
  - `*_index*.bt2` - Bowtie2 index files (after alignment)
  - `*_paired.bam`, `*_singleton.bam` - Intermediate BAMs (after merging)
  - `*_paired.log`, `*_singleton.log` - Log files
- **Keeps**: Final merged BAM file
- **Disk savings**: ~3-5 GB per sample
- **Timing**: Progressive cleanup as each step completes

#### SEMIBIN (Binning)
- **Location**: `modules/local/semibin/main.nf`
- **Removes**:
  - `output/` - Temporary output directory
  - `contig_output/` - Intermediate feature data
  - `*_semibin_bins/` - Temporary bins directory
- **Keeps**: Final bins in `*_semibin_output_bins/`
- **Disk savings**: ~500 MB - 1 GB per sample
- **Timing**: After moving final bins
- **Note**: Already implemented, verified comprehensive

### Phase 3: Nextflow Configuration

#### Cache Strategy
- **File**: `conf/base.config`
- **Setting**: `cache = 'lenient'`
- **Purpose**: Allows resume even when work directories are cleaned
- **Behavior**: Nextflow looks for outputs in `storeDir`/`publishDir` locations instead of work directory

#### Cleanup Profile
- **File**: `nextflow.config`
- **Profile**: `low_disk`
- **Settings**:
  - `cleanup = true` - Enable automatic work directory cleanup
  - `enable_work_cleanup = true` - Pipeline parameter
  - `store_clean_reads = true` - Enable BOWTIE2 storeDir
  - `store_filtered_contigs = true` - Enable BBMAP storeDir
  - `store_refined_bins = true` - Enable METAWRAP storeDir
  - `process.cache = 'lenient'` - Support resume with cleanup

## Usage

### Option 1: Use the low_disk Profile (Recommended)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker,low_disk \
    -resume
```

This automatically enables:
- ✅ Work directory cleanup
- ✅ storeDir for all large intermediates
- ✅ Lenient cache for resume support

### Option 2: Enable via Command Line Parameters

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --enable_work_cleanup \
    --store_clean_reads \
    --store_filtered_contigs \
    --store_refined_bins \
    -profile docker \
    -resume
```

### Option 3: Selective Cleanup

Enable only specific cleanup features:

```bash
# Only clean reads storage
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --store_clean_reads \
    -profile docker

# Only filtered contigs storage
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --store_filtered_contigs \
    -profile docker
```

## Monitoring Disk Usage

A monitoring script is provided to track disk usage during pipeline execution:

```bash
# Start monitoring in background
./bin/monitor_disk_usage.sh ./work ./results 60 &
MONITOR_PID=$!

# Run pipeline
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker,low_disk \
    -resume

# Stop monitoring when done
kill $MONITOR_PID
```

The script generates:
- `disk_usage_monitor.log` - Human-readable log with timestamps
- `disk_usage_monitor.csv` - CSV data for analysis/plotting

### Monitor Script Options

```bash
./bin/monitor_disk_usage.sh [work_dir] [output_dir] [interval_seconds]

# Examples:
./bin/monitor_disk_usage.sh                    # Default: ./work, ./results, 60s
./bin/monitor_disk_usage.sh ./work ./results   # Custom dirs, 60s interval
./bin/monitor_disk_usage.sh ./work ./results 30  # 30s interval
```

## Expected Disk Space Impact

### Per Sample (Typical Metagenomic Dataset)

| Stage | Before | After | Savings |
|-------|--------|-------|---------|
| QC (clean reads) | 10 GB | 2 GB | 8 GB |
| Assembly (contigs) | 8 GB | 3 GB | 5 GB |
| Alignment (BAMs) | 12 GB | 7 GB | 5 GB |
| Binning (bins) | 5 GB | 2 GB | 3 GB |
| **Total per sample** | **35 GB** | **14 GB** | **21 GB** |

### For 10 Samples

| Metric | Without Cleanup | With Cleanup | Reduction |
|--------|----------------|--------------|-----------|
| Peak work directory | 250-300 GB | 80-120 GB | **60-70%** |
| Final output directory | 50 GB | 50 GB | 0% (same) |
| Total disk required | 300-350 GB | 130-170 GB | **50-60%** |

## Timeline Example (3 Samples)

```
Time | Stage              | Action                    | Work Dir | Change
-----|-------------------|---------------------------|----------|--------
T0   | Start             | -                         | 0 GB     | -
T1   | QC running        | FASTP + BOWTIE2           | 60 GB    | +60 GB
T2   | QC complete       | Store clean reads         | 10 GB    | -50 GB ✅
T3   | Assembly running  | MEGAHIT                   | 40 GB    | +30 GB
T4   | Assembly cleanup  | Remove k-mer dirs         | 25 GB    | -15 GB ✅
T5   | BBMAP complete    | Store filtered contigs    | 10 GB    | -15 GB ✅
T6   | Alignment running | BOWTIE2_SAMTOOLS          | 50 GB    | +40 GB
T7   | Alignment cleanup | Remove indices/temp BAMs  | 30 GB    | -20 GB ✅
T8   | Binning running   | MetaBAT2/SemiBin/COMEBin  | 45 GB    | +15 GB
T9   | Binning complete  | Store refined bins        | 10 GB    | -35 GB ✅
T10  | Pipeline complete | Final cleanup             | 5 GB     | -5 GB ✅
```

**Key Point**: Disk space is freed progressively throughout execution, preventing the pipeline from running out of space.

## Resume Functionality

The optimization preserves Nextflow's resume capability:

### How Resume Works with Cleanup

1. **storeDir outputs** are permanent and never deleted
2. **Cache strategy** is set to `lenient`
3. **Nextflow checks** for outputs in storeDir locations
4. **If outputs exist**, task is cached (not re-run)
5. **If outputs missing**, task is re-executed

### Testing Resume

```bash
# Run 1: Start pipeline
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker,low_disk

# Interrupt it (Ctrl+C) during execution

# Run 2: Resume from where it stopped
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker,low_disk \
    -resume

# Check cached processes
grep "Cached" .nextflow.log | wc -l
```

### What Gets Cached

- ✅ Processes with storeDir outputs (BOWTIE2, BBMAP, METAWRAP)
- ✅ Processes with publishDir outputs (all other processes)
- ✅ Processes whose work dirs were cleaned but outputs exist

### What Gets Re-run

- ❌ Processes whose outputs were manually deleted
- ❌ Processes with changed inputs or parameters
- ❌ Processes with modified scripts

## Troubleshooting

### Issue: Pipeline runs out of disk space

**Cause**: Cleanup not enabled or insufficient disk space for peak usage

**Solution**:
1. Ensure you're using `-profile low_disk`
2. Check that storeDir parameters are enabled
3. Monitor disk usage with the monitoring script
4. Consider reducing number of parallel samples

### Issue: Resume not working after cleanup

**Cause**: Outputs not found in expected locations

**Solution**:
1. Verify `cache = 'lenient'` is set in `conf/base.config`
2. Check that storeDir outputs exist in output directory
3. Ensure output directory path hasn't changed
4. Check `.nextflow.log` for cache lookup messages

### Issue: Work directory still growing

**Cause**: Cleanup happens after process completion, not during

**Solution**:
1. This is expected - work dir grows during process execution
2. Cleanup occurs when process completes successfully
3. Monitor with the disk usage script to see cleanup patterns
4. Consider reducing parallel process count with `max_cpus`

### Issue: Outputs missing from results directory

**Cause**: storeDir only stores specific intermediate files

**Solution**:
1. Check that you're looking in the correct subdirectory:
   - Clean reads: `output/clean_reads/{sample_id}/`
   - Filtered contigs: `output/assembly/{sample_id}/`
   - Refined bins: `output/bins/{sample_id}/`
2. Other outputs use standard publishDir locations
3. Check process-specific output directories

## Configuration Parameters

### Pipeline Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `enable_work_cleanup` | `false` | Enable automatic work directory cleanup |
| `store_clean_reads` | `false` | Use storeDir for clean reads (BOWTIE2) |
| `store_filtered_contigs` | `false` | Use storeDir for filtered contigs (BBMAP) |
| `store_refined_bins` | `false` | Use storeDir for refined bins (METAWRAP) |

### Process Settings

| Setting | Value | Description |
|---------|-------|-------------|
| `cache` | `'lenient'` | Cache strategy for resume support |
| `cleanup` | `true` (in low_disk profile) | Enable work directory cleanup |

## Best Practices

1. **Always use `-resume`** when restarting failed runs
2. **Monitor disk usage** during first run to understand patterns
3. **Use `low_disk` profile** on systems with limited disk space
4. **Keep output directory** on a filesystem with sufficient space
5. **Don't manually delete** work directories during execution
6. **Check logs** if resume doesn't work as expected
7. **Test with small dataset** first to verify cleanup behavior

## Technical Details

### publishDir Modes

| Mode | Behavior | Disk Impact | Use Case |
|------|----------|-------------|----------|
| **copy** | Copy files to output | No space freed | Default, safe for all outputs |
| **move** | Move files to output | Immediate space freed | Large intermediates, cleanup enabled |
| **symlink** | Create symbolic links | Minimal space | Read-only access |
| **rellink** | Create relative links | Minimal space | Portable links |

### Why mode 'move' Works for Cleanup

- **Immediate**: Files are moved as soon as the process completes
- **No duplication**: Files exist only in output directory, not in work directory
- **Resume compatible**: Nextflow finds outputs in publishDir location with lenient cache
- **Selective**: Only enabled when cleanup parameters are true

### Cache Strategies

| Strategy | Behavior | Use Case |
|----------|----------|----------|
| `standard` | Check work dir only | Default, no cleanup |
| `lenient` | Check work dir + publishDir/storeDir | With cleanup enabled |
| `deep` | Check all inputs deeply | Strict reproducibility |

### Cleanup Timing

- **Internal cleanup**: During process execution (rm commands in script)
- **storeDir cleanup**: Immediately after outputs are moved
- **Nextflow cleanup**: After downstream processes no longer need the work dir

## Version History

- **v1.0.0** - Initial disk optimization implementation
  - Phase 1: storeDir for BOWTIE2, BBMAP, METAWRAP
  - Phase 2: Internal cleanup for MEGAHIT, BOWTIE2_SAMTOOLS
  - Phase 3: Nextflow configuration and monitoring script

## Support

For issues or questions about disk space optimization:
1. Check this documentation first
2. Review `.nextflow.log` for detailed execution information
3. Use the monitoring script to track disk usage patterns
4. Open an issue on the GitHub repository with logs and disk usage data
