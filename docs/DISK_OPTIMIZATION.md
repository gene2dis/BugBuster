# BugBuster Pipeline - Disk Space

This document describes what BugBuster does to limit disk use, what the `low_disk` profile
does (and does not do), and how to run on a machine with limited disk space.

## Overview

Disk use comes from three places, which behave differently:

1. **The work directory** (`work/`, or `-work-dir`): every task's inputs, intermediates and
   outputs. This is where peak usage happens, and it is also the `-resume` cache.
2. **The results directory** (`--output`): copies of the published outputs.
3. **The databases** (`--databases_dir`, default `<output>/../databases`): reference
   databases, kept across runs.

BugBuster reduces disk use in two ways:

1. **In-task cleanup**, always on: a few large processes delete their own temporary files
   before they finish (see below). This lowers the size of their task directories in every
   run, with or without `low_disk`.
2. **Work-dir deletion after a successful run**: Nextflow's `cleanup = true`, set by the
   `low_disk` profile. It runs **once, when the run completes successfully**. It does not
   free space while the run is in progress, so it does **not** lower peak disk use.

All publishing uses standard `publishDir` with `mode: 'copy'` (`publish_dir_mode`).
Publishing never removes files from the work directory.

## In-Task Cleanup (all runs)

These processes delete temporary files inside their own task directory before finishing:

#### MEGAHIT (Assembly)
- **Location**: `modules/local/megahit/main.nf`
- **Removes**: `intermediate_contigs/`, `kmer_k*/`, `*.tmp`, `checkpoints.txt`
- **Keeps**: final contigs and log files

#### BOWTIE2_SAMTOOLS (Alignment)
- **Location**: `modules/local/bowtie2_samtools/main.nf`
- **Removes**: the per-task Bowtie2 index (`*_index*.bt2`) after alignment, the
  intermediate `*_paired.bam` / `*_singleton.bam` after merging, and their logs
- **Keeps**: the final merged BAM

#### SEMIBIN (Binning)
- **Location**: `modules/local/semibin/main.nf`
- **Removes**: `output/`, `contig_output/` and the temporary `*_semibin_bins/` directory
- **Keeps**: the final bins in `*_semibin_output_bins/`

Because these files are removed inside the task, the task's outputs (and therefore
`-resume`) are unaffected.

## The `low_disk` Profile

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker,low_disk
```

It sets:

- `cleanup = true`: when the run **completes successfully**, Nextflow deletes the work
  directories of the tasks that run executed.
- `store_clean_reads = true`: the decontaminated reads are published to
  `<output>/clean_reads/{sample_id}/`, so they survive the work-dir deletion.
- `process.cache = 'lenient'`: the same cache mode every run already uses (from
  `conf/base.config`).

What this means in practice:

- **Peak disk use is unchanged.** The work directory grows during the run exactly as it
  does without the profile. Size the work-dir filesystem for the full run.
- **A completed `low_disk` run cannot be resumed.** Its work directory is gone, so
  re-running with a changed or added option recomputes everything.
- **An interrupted or failed `low_disk` run can be resumed.** Cleanup only runs on
  success, so the work directory is still there.
- **After a resumed run succeeds, some task directories are left behind.** Cleanup deletes
  only the tasks the final run executed. Tasks reused from the cache, and tasks aborted by
  the interruption, stay in `work/`. Delete `work/` once you no longer need to resume.
  `nextflow clean -f` is not enough here: it removes only the last run's tasks, and a
  task directory created just before the interruption may not be recorded in any run.
- The final results directory is the same with or without the profile, plus
  `clean_reads/`.

> **Note**: work-dir cleanup is a Nextflow config setting (`cleanup = true`), not a
> pipeline parameter, so there is no `--flag` for it. Use `-profile low_disk`, or add
> `cleanup = true` to a custom config passed with `-c`.

> **Note: databases are unaffected.** Reference databases live at
> `--databases_dir` (default `<output>/../databases`), outside the work dir,
> and `low_disk` does nothing about them. The functional annotation branches in
> particular add large permanent databases (eggNOG 7 ~44 GB, dbCAN ~7.4 GB,
> Bakta full ~31.9 GB / light ~1.3 GB download, WoLr2 ~94 GB for the
> read-level Woltka branch, SUPER-FOCUS DB_90 ~0.7-0.9 GB download for the
> read-level SUPER-FOCUS branch, HUMAnN ~71 GB with the full ChocoPhlAn or
> ~33 GB EC-filtered for the read-level HUMAnN branch). They also pair badly with
> `low_disk`: a completed run cannot be resumed, so changing an option afterwards
> repeats the long eggNOG-mapper and read-level steps (the pipeline warns about this
> combination at launch). The Bakta light database is never selected automatically:
> `--bakta_db v6.0-light` is always an explicit choice, recorded in provenance.

To publish the clean reads without deleting the work directory, use `--store_clean_reads`
on its own:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --store_clean_reads \
    -profile docker \
    -resume
```

## Running With Limited Disk Space

Peak usage is in the work directory, so:

1. **Put the work directory on the largest filesystem available**:
   `-work-dir /scratch/bugbuster_work`. The results and databases directories can live
   elsewhere.
2. **Run fewer tasks at once.** Fewer concurrent tasks means fewer task directories at
   their largest at the same time. With the local executor, lower `--max_cpus` (it is also
   the executor's CPU pool); on a cluster, lower the executor's `queueSize` in a `-c` config.
3. **Split large studies** into batches of samples (each batch its own run and work dir).
   Co-assembly needs all samples in one run.
4. **Free space between runs**: once a run's results are final, delete its work directory.
   This is what `low_disk` does automatically on success. (`nextflow clean -f` removes
   only the last run's tasks, so after resumed runs it leaves the earlier ones behind.)
5. **Monitor usage** with the monitoring script below to see where the peak is.

## Monitoring Disk Usage

A monitoring script tracks disk usage during pipeline execution:

```bash
# Start monitoring in background
./bin/monitor_disk_usage.sh ./work ./results 60 &
MONITOR_PID=$!

# Run pipeline
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker \
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

## Resume Functionality

### How Resume Interacts with Cleanup

1. `-resume` relies on the **work directory**: a task is cached only while its work dir
   (and outputs) still exist.
2. The cache mode is `lenient` in every run (`conf/base.config`): input files are matched
   by path and size, ignoring timestamps.
3. `cleanup = true` (the `low_disk` profile) deletes the work directory **at the end of a
   successful run**, so a completed `low_disk` run cannot be resumed; an interrupted one
   can.

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

# Check cached processes: tasks reused from the first run show status CACHED
nextflow log <run name> -f name,status    # run names: nextflow log
```

### What Gets Cached

- ✅ Processes whose work dirs still exist (normal `-resume` behavior)

> **Note**: publishing is not caching. Published outputs (clean reads,
> filtered contigs, refined bins) are copies for the user; every task
> re-runs whenever its inputs or parameters change, regardless of what has
> been published. `-resume` (the work dir) is the cache.

### What Gets Re-run

- ❌ Processes whose work dirs were deleted (including by `low_disk` after a completed run)
- ❌ Processes with changed inputs or parameters
- ❌ Processes with modified scripts

## Troubleshooting

### Issue: Pipeline runs out of disk space

**Cause**: The work-dir filesystem is smaller than the run's peak usage. `low_disk` does
not help here: it deletes the work directory only after the run succeeds.

**Solution**:
1. Start the run again with `-work-dir` on a larger filesystem. Do not count on resuming
   from a work directory copied elsewhere: task directories refer to each other by
   absolute path
2. Run fewer tasks at once, or split the samples into batches (see "Running With Limited
   Disk Space")
3. Monitor disk usage with the monitoring script

### Issue: Resume not working after cleanup

**Cause**: The work directory was deleted. A completed `low_disk` run deletes it.

**Solution**:
1. Verify the `work/` directory still exists
2. If it is gone, the run has to be recomputed. Use `low_disk` only for runs you will not
   need to resume or extend.

### Issue: Work directory not empty after a successful `low_disk` run

**Cause**: The run was resumed. Cleanup deletes only the tasks the final run executed;
tasks reused from the cache and tasks aborted by the earlier interruption remain.

**Solution**: delete `work/` once you no longer need to resume. `nextflow clean -f` alone
removes only the last run's tasks; `nextflow clean -f <run name>` (names from
`nextflow log`) removes an earlier run's recorded tasks, but a task directory created just
before an interruption may not be recorded in any run.

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
| `cache` | `'lenient'` (all runs, `conf/base.config`) | Match input files by path and size, ignoring timestamps |
| `cleanup` | `true` (in the `low_disk` profile) | Delete the work directory after a successful run |

## Best Practices

1. **Always use `-resume`** when restarting failed runs
2. **Size the work-dir filesystem for the peak**, and put it on the largest disk available
3. **Monitor disk usage** during a first run to learn where the peak is
4. **Use `low_disk` for final runs** you will not need to resume or extend
5. **Keep the output directory** on a filesystem with sufficient space
6. **Don't manually delete** work directories during execution

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
| `standard` | Match input files by path, size and last-modified time | Nextflow default |
| `lenient` | Match input files by path and size only (tolerates touched files) | What BugBuster uses |
| `deep` | Match input files by content | Strict reproducibility |

### Cleanup Timing

- **In-task cleanup**: during process execution (rm commands in the process script), every run
- **Nextflow cleanup** (`cleanup = true`, `low_disk` only): once, at the end of a successful run

## Version History

- **v2.0.0dev** - Current: this document corrected to the measured behaviour (cleanup does not lower peak usage)
- **v1.1** - `low_disk` profile (work-dir deletion after success + clean-reads publishing); storeDir replaced by conditional publishDir
- **v1.0.0** - Initial disk optimization implementation
  - Phase 1: storeDir for BOWTIE2, BBMAP, METAWRAP
  - Phase 2: Internal cleanup for MEGAHIT, BOWTIE2_SAMTOOLS
  - Phase 3: Nextflow configuration and monitoring script

## Support

For issues or questions about disk space:
1. Check this documentation first
2. Review `.nextflow.log` for detailed execution information
3. Use the monitoring script to track disk usage patterns
4. Open an issue on the GitHub repository with logs and disk usage data
