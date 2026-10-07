# BugBuster Pipeline - Disk Space Optimization

This document describes the disk space optimization features implemented in BugBuster to prevent running out of disk space during pipeline execution.

## Overview

The pipeline implements a **progressive cleanup strategy** that frees disk space **during execution**, not after. This is achieved through two complementary mechanisms:

1. **Internal process cleanup** - Remove temporary files within processes as soon as they're no longer needed
2. **Nextflow cleanup configuration** - Enable automatic work directory cleanup

All publishing uses standard `publishDir` with `mode: 'copy'` (`publish_dir_mode`). Publishing never removes files from the work directory; disk is reclaimed by the internal cleanup steps below and by Nextflow's `cleanup` setting.

## Implementation Details

### Phase 1: Key Output Publishing

- **BOWTIE2_DECONTAMINATE (clean reads)**: published to `output/clean_reads/{sample_id}/` only when `--store_clean_reads` is set (enabled by the `low_disk` profile). Configured in `config/modules.config`.
- **BBMAP (filtered contigs)**: always published to `output/03_assembly/per_sample/{sample_id}/` (or `03_assembly/coassembly/`) — `*_filtered_contigs.fa` and `*_contig.stats`.
- **METAWRAP (refined bins)**: always published to `output/04_binning/per_sample/{sample_id}/refined_bins/` (or `04_binning/coassembly/refined_bins/`) — the `*_metawrap_*_bins/` directory.

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
  - `store_clean_reads = true` - Publish clean reads to `<output>/clean_reads/` (publishDir)
  - `process.cache = 'lenient'` - Support resume with cleanup

## Usage

### Use the low_disk Profile

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker,low_disk
```

This automatically enables:
- ✅ Work directory cleanup (`cleanup = true` — the work dir is deleted after a successful run, so **low_disk runs are not resumable**)
- ✅ Clean-reads publishing (`store_clean_reads`)
- ✅ Lenient cache

> **Note**: Work-dir cleanup is a Nextflow config setting (`cleanup = true`), not a pipeline parameter — there is no `--flag` for it. Use `-profile low_disk`, or add `cleanup = true` to a custom config passed with `-c`.

> **Note — databases are unaffected**: reference databases live at
> `--databases_dir` (default `<output>/../databases`), outside the work dir,
> and `low_disk` does nothing about them. The functional annotation branches in
> particular add large permanent databases (eggNOG 7 ~44 GB, dbCAN ~7.4 GB,
> Bakta full ~31.9 GB / light ~1.3 GB download, WoLr2 ~94 GB for the
> read-level Woltka branch, SUPER-FOCUS DB_90 ~0.7-0.9 GB download for the
> read-level SUPER-FOCUS branch) —
> and pair badly with `low_disk`, since the long eggNOG runs are exactly
> where `-resume` matters most (the pipeline warns about this combination at
> launch). The Bakta light database is never selected automatically:
> `--bakta_db v6.0-light` is always an explicit choice, recorded in provenance.

To publish clean reads without the cleanup trade-off, use `--store_clean_reads` on its own:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --store_clean_reads \
    -profile docker \
    -resume
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

## Resume Functionality

### How Resume Interacts with Cleanup

1. `-resume` relies on the **work directory**: a task is cached only while its work dir (and outputs) still exist
2. **Cache strategy** is set to `lenient` to tolerate internal cleanup of temporary files inside completed task dirs
3. `cleanup = true` (the `low_disk` profile) deletes the work directory **at the end of a successful run**, so a completed `low_disk` run cannot be resumed; an interrupted one can

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

- ✅ Processes whose work dirs still exist (normal `-resume` behavior)

> **Note**: publishing is not caching. Published outputs (clean reads,
> filtered contigs, refined bins) are copies for the user; every task
> re-runs whenever its inputs or parameters change, regardless of what has
> been published. `-resume` (the work dir) is the cache.

### What Gets Re-run

- ❌ Processes whose outputs were manually deleted
- ❌ Processes with changed inputs or parameters
- ❌ Processes with modified scripts

## Troubleshooting

### Issue: Pipeline runs out of disk space

**Cause**: Cleanup not enabled or insufficient disk space for peak usage

**Solution**:
1. Ensure you're using `-profile low_disk`
2. Monitor disk usage with the monitoring script
3. Consider reducing number of parallel samples

### Issue: Resume not working after cleanup

**Cause**: Outputs not found in expected locations

**Solution**:
1. Verify the `work/` directory still exists (it is deleted at the end of a completed `low_disk` run)
2. Verify `cache = 'lenient'` is set in `conf/base.config`
3. Check `.nextflow.log` for cache lookup messages

### Issue: Work directory still growing

**Cause**: Cleanup happens after process completion, not during

**Solution**:
1. This is expected - work dir grows during process execution
2. Cleanup occurs when process completes successfully
3. Monitor with the disk usage script to see cleanup patterns
4. Consider reducing parallel process count with `max_cpus`

### Issue: Outputs missing from results directory

**Solution**:
1. Check that you're looking in the correct subdirectory:
   - Clean reads: `output/clean_reads/{sample_id}/` (only with `--store_clean_reads`)
   - Filtered contigs: `output/03_assembly/per_sample/{sample_id}/` (or `03_assembly/coassembly/`)
   - Refined bins: `output/04_binning/per_sample/{sample_id}/refined_bins/` (or `04_binning/coassembly/refined_bins/`)
2. Other outputs use their numbered publishDir locations (see `config/modules.config`)

## Configuration Parameters

### Pipeline Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `store_clean_reads` | `false` | Publish clean reads to the output dir (BOWTIE2_DECONTAMINATE, publishDir) |

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
| **copy** | Copy files to output | No space freed | Default; what BugBuster uses (`publish_dir_mode`) |
| **move** | Move files to output | Immediate space freed | NOT used: it removes outputs downstream tasks still stage from the work dir and breaks `-resume` |
| **symlink** | Create symbolic links | Minimal space | Read-only access |
| **rellink** | Create relative links | Minimal space | Portable links |

### Cache Strategies

| Strategy | Behavior | Use Case |
|----------|----------|----------|
| `standard` | Check work dir only | Default, no cleanup |
| `lenient` | Match inputs by path/size only (tolerates touched files) | With internal cleanup enabled |
| `deep` | Check all inputs deeply | Strict reproducibility |

### Cleanup Timing

- **Internal cleanup**: During process execution (rm commands in script)
- **Nextflow cleanup** (`cleanup = true`): at the end of a successful run

## Version History

- **v1.1.0dev** - Current: `low_disk` profile (work-dir cleanup + clean-reads publishing); storeDir replaced by conditional publishDir
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
