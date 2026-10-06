# Troubleshooting Guide

This guide covers common issues and solutions when running BugBuster.

## Table of Contents

- [Installation Issues](#installation-issues)
- [Input/Samplesheet Errors](#inputsamplesheet-errors)
- [Resource Errors](#resource-errors)
- [Container Issues](#container-issues)
- [Database Issues](#database-issues)
- [Decontamination Issues](#decontamination-issues)
- [Cloud Execution Issues](#cloud-execution-issues)
- [Output Issues](#output-issues)
- [Debug Mode](#debug-mode)
- [Quick Fixes Checklist](#quick-fixes-checklist)
- [Getting Help](#getting-help)

---

## Installation Issues

### Nextflow version too old

**Error:**
```
ERROR: Nextflow version 25.10.0 or later is required
```

**Solution:**
```bash
# Update Nextflow
nextflow self-update

# Or install specific version
curl -s https://get.nextflow.io | bash
./nextflow self-update
```

### Java not found

**Error:**
```
ERROR: Cannot find Java or it's a wrong version
```

**Solution:**
```bash
# Install Java 11 or later
sudo apt-get install openjdk-11-jdk

# Or use SDKMAN
curl -s "https://get.sdkman.io" | bash
sdk install java 11.0.21-tem
```

---

## Input/Samplesheet Errors

### Invalid samplesheet format

**Error:**
```
ERROR: Invalid samplesheet: 'sample' column is empty
```

**Solution:**
1. Ensure CSV has header: `sample,r1,r2,s`
2. Check for extra spaces or special characters
3. Use absolute paths for files
4. Validate with:
   ```bash
   head -5 samplesheet.csv
   ```

### Files not found

**Error:**
```
ERROR: Cannot find file: /path/to/reads.fastq.gz
```

**Solution:**
1. Use absolute paths in samplesheet
2. Check file permissions: `ls -la /path/to/reads.fastq.gz`
3. Ensure files are not compressed with unsupported format

### Invalid sample name

**Error:**
```
Invalid sample name '...': use only letters, digits, underscore, dot or hyphen
```

**Solution:**
- Sample names become file and directory names throughout the pipeline
- Use only `A-Za-z0-9`, `_`, `.` or `-`, starting with a letter or digit
- Replace spaces with underscores

### Unrecognised parameter

**Error:**
```
* --<name>: expected type: ... (unrecognised parameter)
```
or a validation error naming a parameter you passed.

**Solution:**
- Parameters are validated against `nextflow_schema.json` at startup; any parameter the pipeline does not declare aborts the run
- Check the spelling against [`parameters.md`](parameters.md)
- If you are following instructions written for an older release, the parameter may have been removed or renamed (e.g. `--kraken_db_used`, `--sourmash_db_name`, `--validationShowHiddenParams`, `--enable_work_cleanup`, and all `--mmseqs_*` parameters no longer exist; `--bbmap_lenght` is now `--bbmap_length`; `--custom_phiX_index` and `--custom_bowtie_host_index` were replaced by `--custom_decontamination_index` / `--custom_phiX_fasta` / `--custom_host_fasta`)

---

## Resource Errors

### Out of memory

**Error:**
```
Process exceeded memory limit
exitCode: 137
```

**Solution:**
1. Increase memory limits:
   ```bash
   nextflow run main.nf --max_memory '256.GB' ...
   ```

2. For specific processes, edit `conf/base.config`:
   ```groovy
   process {
       withName: 'MEGAHIT.*' {
           memory = '128.GB'
       }
   }
   ```

3. Use a larger instance type (cloud) or node (HPC)

### Out of disk space

**Error:**
```
No space left on device
```

**Solution:**
1. Clean work directory:
   ```bash
   nextflow clean -f
   ```

2. Use external storage for work directory:
   ```bash
   nextflow run main.nf -work-dir /scratch/work ...
   ```

3. Enable cleanup on success in config:
   ```groovy
   cleanup = true
   ```

### Process timeout

**Error:**
```
Process exceeded time limit
```

**Solution:**
1. Increase time limits:
   ```bash
   nextflow run main.nf --max_time '480.h' ...
   ```

2. For specific processes:
   ```groovy
   process {
       withName: 'GTDB_TK.*' {
           time = '168.h'
       }
   }
   ```

---

## Container Issues

### Docker permission denied

**Error:**
```
permission denied while trying to connect to the Docker daemon
```

**Solution:**
```bash
# Add user to docker group
sudo usermod -aG docker $USER

# Log out and back in, or run:
newgrp docker
```

### Singularity image pull failed

**Error:**
```
FATAL: Unable to pull docker://quay.io/biocontainers/...
```

**Solution:**
1. Check network connectivity
2. Set cache directory:
   ```bash
   export SINGULARITY_CACHEDIR=/path/to/cache
   export NXF_SINGULARITY_CACHEDIR=/path/to/cache
   ```

3. Increase pull timeout in config:
   ```groovy
   singularity {
       pullTimeout = '60 min'
   }
   ```

### Container not found

**Error:**
```
Unable to find image 'container:tag' locally
```

**Solution:**
1. Check internet connectivity
2. Verify container exists (look up the exact tag in the module's `main.nf`):
   ```bash
   docker pull quay.io/biocontainers/fastp:<tag>
   ```

3. Use alternative registry if blocked

### Functional annotation aborts under docker/podman

**Error:**
```
--contig_level_functional requires a singularity or apptainer container engine ...
```

**Cause:** eggNOG-mapper v3 is in beta and its authors publish only an Apptainer
image — there is no docker image or biocontainer yet. The pipeline rejects real
runs of the functional annotation branch under docker/podman at launch instead
of failing mid-run on a container pull.

**Solution:** run with a singularity or apptainer profile:
```bash
nextflow run main.nf ... --contig_level_functional -profile apptainer
```
Nextflow uses one container engine per run, so this switches the whole run —
not just the eggNOG steps — to singularity/apptainer. Every other tool runs
from the same pinned images (automatically converted on first use), so results
are identical; only the first-run conversion time differs. The docker-based
cloud profiles (`aws`, `gcp`, `azure`) cannot be combined with
`--contig_level_functional` during the beta. Stub runs (`-stub`) are
unaffected. Docker and cloud support return once eggNOG-mapper v3.0.0 final is
released on bioconda. This guard applies only to the contig branch
(`--contig_level_functional`): the MAG branch (`--mag_level_functional`,
Bakta) uses a normal biocontainer and runs on any engine.

> Related: the upstream eggNOG-mapper image ships without `procps`, which would
> normally abort Nextflow tasks with `Command 'ps' required by nextflow to
> collect task metrics cannot be found`. The pipeline works around this with a
> delegating `ps` shim in `bin/` — do not remove `bin/ps` while the functional
> branch uses the upstream image.

### AGGREGATE_FUNCTIONS fails with "annotations column header does not match the verified ... layout" or "not among the layouts this parser was verified against"

**Error:**
```
Error: ...emapper.annotations: annotations column header does not match the verified emapper-3.0.0-beta6 layout ...
Error: eggnog-mapper version '<X>' is not among the layouts this parser was verified against ...
```

**Cause:** this is a deliberate guard, not a bug. eggNOG-mapper v3 is a beta and
its output schema can change between releases, so the functional aggregation
step (`bin/aggregate_functions.py`) only parses annotations from eggNOG-mapper
versions whose column layout it was explicitly verified against, and checks the
file's column header against that layout. Silently misparsing a shifted column
into the wrong ontology would corrupt every downstream table, which is exactly
what this refuses to do.

**Solution:** this error appears when the pinned eggNOG-mapper version and the
aggregation script get out of sync (for example after re-pinning the eggNOG
modules to a newer release without updating the parser). Re-verify the new
version's column layout on real output, then update `KNOWN_EMAPPER_VERSIONS`
and, if the columns changed, `EXPECTED_ANNOTATION_COLUMNS` in
`bin/aggregate_functions.py` — and extend
`tests/bin/test_aggregate_functions.sh` for the new layout. Do not bypass the
check by editing the annotations file.

The same guard exists for the run_dbcan overview (`overview column header does
not match the verified run_dbcan v5 layout ...` / `run_dbcan version '<X>' is
not among the layouts ...`): after re-pinning the `RUN_DBCAN` module, re-verify
the overview columns on real output and update `KNOWN_DBCAN_VERSIONS` and, if
needed, `DBCAN_OVERVIEW_COLUMNS` in `bin/aggregate_functions.py`, plus the
harness tests.

---

### `--mag_level_functional` rejected at launch

**Error:**
```
--mag_level_functional requires --include_binning (Bakta annotates the refined bins)
--mag_level_functional requires at least two --binners: MetaWRAP refinement and its completeness/contamination quality filter only run with >= 2 binners ...
```

**Cause:** deliberate validation, not a bug. Bakta annotates the refined bins,
so binning must be enabled — and the pipeline's only bin quality filter is
MetaWRAP refinement (`-c`/`-x`, defaults 50/10), which only runs when at least
two binners are selected. With a single binner the bins channel carries raw,
unfiltered binner output, and the MAG branch refuses to annotate unfiltered
bins.

**Solution:** enable binning with two or more binners:
```bash
nextflow run main.nf ... --include_binning --binners semibin,metabat2 --mag_level_functional
```
Note that Bakta annotates every MetaWRAP-refined bin; it does not further
filter on CheckM2 scores, and it does not detect or exclude
archaeal/eukaryotic/viral bins — check the GTDB-Tk bin taxonomy report before
interpreting annotations of non-bacterial bins.

### FORMAT_BAKTA_DB fails: schema major / amrfinderplus-db errors

**Error:**
```
ERROR: bakta_db/version.json reports schema major '<X>', expected 6 (required by Bakta 1.12.x)
ERROR: bakta_db/amrfinderplus-db missing or empty after extraction (corrupt or truncated download?)
```

**Cause:** the Bakta database download is hard-verified after extraction.
Bakta 1.12.x only accepts database schema 6 (the pinned Zenodo v6.0 release),
and the release tarball bundles its AMRFinderPlus database — a missing or
empty `amrfinderplus-db/` means a corrupt or truncated download, which must
fail rather than be silently "repaired".

**Solution:** for the schema error, check that `--custom_bakta_db` points at
an extracted v6.0 (schema 6) database, not an older release. For the
amrfinderplus-db error, delete the partial `databases/bakta/` directory and
re-run — the download resumes from the pinned Zenodo URL. The same checks
apply to a custom database directory.

### BAKTA_BAKTA fails mid-annotation with "amrfinder error! error code: 1"

**Error (after ~10 minutes of successful annotation, at the AMR expert step):**
```
Exception: amrfinder error! error code: 1. Please, try 'amrfinder_update --force_update --database .../amrfinderplus-db' ...
```
and the per-bin `<sample>_<bin>.log` shows
`Software requires database version at least <date>` when amrfinder is run by hand.

**Cause:** the official Bakta DB v6.0 tarball bundles a 2024-era AMRFinderPlus
database, but the AMRFinderPlus **program** inside the pinned bakta container
is newer and refuses databases older than its minimum. This bites
`--custom_bakta_db` directories extracted straight from the official tarball.
The pipeline's own download path (`FORMAT_BAKTA_DB`) refreshes the
AMRFinderPlus component during provisioning, so auto-downloaded databases are
not affected.

**Solution:** refresh the AMRFinderPlus component inside your custom database
once (~300 MB from NCBI; adds a new dated subdirectory and repoints the
`latest` symlink — the rest of the database is untouched):
```bash
docker run --rm -u $(id -u):$(id -g) \
    -v /path/to/bakta_db/amrfinderplus-db:/amrdb \
    quay.io/biocontainers/bakta:1.12.1--pyhdfd78af_0 \
    amrfinder_update --force_update --database /amrdb
```
Then re-run with `-resume` (the Bakta tasks re-run because their database
input changed; everything upstream stays cached).

### MICROBECENSUS task fails / `cpge` columns are empty / `ags_and_ge.tsv` says `unavailable`

**Symptom:** the run completes, but a `MICROBECENSUS` task shows as failed
(ignored) in the log, and in `07_functional_annotation/summary/` the affected
sample has empty `cpge` fields, a blank column in the
`function_wide_*_cpge.tsv` matrices, and `status = unavailable` in
`ags_and_ge.tsv`.

**Cause:** this is the designed behavior, not a crash. MicrobeCensus failure
is deliberately non-fatal: the sample falls back to TPM-only normalization and
the run continues. The common reasons it fails:

1. **Reads shorter than 50 bp** — MicrobeCensus has no model below 50 bp and
   exits with `Cannot compute AGS using reads shorter than 50 bp`.
2. **Too few marker-gene hits** — very small datasets (subsampled tests,
   shallow runs) exit with `No hits to marker proteins - cannot estimate
   genome size`. Estimates need a few hundred thousand reads to be
   meaningful.
3. **Implausible estimate rejected** — the module refuses an average genome
   size outside 0.5–20 Mb (the upstream tool can occasionally emit silent
   garbage estimates; see MicrobeCensus issue #36).

The exact reason is in the failed task's `.command.log` under the work
directory printed in the Nextflow log.

**Solution:** nothing to fix for the run itself — TPM tables are complete and
valid. If you need CPGE for that sample, address the cause (deeper
sequencing, reads ≥ 50 bp) and re-run; `--microbecensus false` turns the step
off entirely.

**Note:** MicrobeCensus is invoked through `bin/run_microbe_census_py3fix.py`.
Both published biocontainers of the unmaintained upstream tool are broken
(the Python-3 image lacks an unreleased upstream str/bytes fix — upstream
issue #32 — and the Python-2 image cannot load its bundled RAPsearch2
binary); the shim patches the one broken function in the pinned Python-3
image and delegates everything else to the stock tool. Do not remove the shim
while the pin is `microbecensus:1.1.1--pyhca03a8a_2`.

---

## Database Issues

### Database download failed

**Error:**
```
ERROR: Failed to download database from URL
```

**Solution:**
1. Check internet connectivity
2. Use pre-downloaded databases:
   ```bash
   nextflow run main.nf \
       --custom_kraken_db /path/to/kraken_db \
       --custom_checkm2_db /path/to/checkm2_db \
       ...
   ```

3. Download manually and specify path

### Database checksum mismatch

**Error:**
```
ERROR: Database checksum verification failed
```

**Solution:**
1. Delete corrupted download and retry
2. Download from alternative source
3. Use custom database path

### GTDB-TK memory error

**Error:**
```
MemoryError in GTDB-TK
```

**Solution:**
1. GTDB-TK requires ~256GB RAM for full database
2. Increase memory allocation:
   ```groovy
   process {
       withName: 'GTDB_TK.*' {
           memory = '256.GB'
       }
   }
   ```

### GTDB-TK rejects the reference data (wrong release)

**Symptom:** `GTDB_TK_BATCH` fails at startup complaining about the reference
data version or missing files under `GTDBTK_DATA_PATH`.

**Cause:** the pipeline pins GTDB-Tk 2.7.2, which accepts **only GTDB R232**
reference data. Older R220/R226 packages — including a `databases/gtdbtk/`
directory cached by a previous pipeline version — require GTDB-Tk ≤ 2.6.1 and
no longer work.

**Solution:** delete the cached `<databases_dir>/gtdbtk/` directory and let the
pipeline download the R232 package (~61 GB), or pass an R232 directory you
already have with `--custom_gtdbtk_db /path/to/release232` (the directory that
directly contains `markers/`, `skani/`, `taxonomy/`, ...).

---

## Decontamination Issues

### Process BOWTIE2_BUILD_COMBINED failed

**Possible causes:** FASTA files not found, insufficient memory, or corrupted FASTA files.

**Solutions:**

1. Verify the phiX/host FASTA files exist and are valid:
   ```bash
   ls -lh /path/to/host.fasta
   zcat -f /path/to/host.fasta | head -n 5
   gunzip -t /path/to/host.fasta.gz   # for gzipped input
   ```
2. Increase memory for the index build:
   ```bash
   nextflow run main.nf ... --max_memory 64.GB
   ```

### Process BOWTIE2_DECONTAMINATE failed

**Possible causes:** index not built correctly, corrupted read files, or insufficient disk space.

**Solutions:**

1. Check the index files (under your databases directory, `<output>/../databases` by default):
   ```bash
   ls -lh <databases_dir>/bowtie_index/contaminants_index/
   # Should list contaminants.1.bt2 ... contaminants.rev.2.bt2
   ```
2. Verify read files decompress cleanly: `zcat sample_R1.fastq.gz | head -n 4`
3. Check disk space: `df -h .`

### Multiple host genomes not being used

**Symptom:** only one genome is used for decontamination.

`custom_host_fasta` takes a **comma-separated string** — not a YAML list, and
with no spaces after the commas:

```yaml
# Correct - comma-separated string
custom_host_fasta: "/path/file1.fasta,/path/file2.fasta,/path/file3.fasta"

# Incorrect - spaces after commas
custom_host_fasta: "/path/file1.fasta, /path/file2.fasta"

# Incorrect - YAML list (not supported)
custom_host_fasta:
  - "/path/file1.fasta"
  - "/path/file2.fasta"
```

### Slow index building or decontamination

- Pre-build the index once and reuse it:
  ```bash
  cat phix.fasta host.fasta > contaminants.fasta
  bowtie2-build --threads 16 contaminants.fasta contaminants_index/contaminants

  nextflow run main.nf ... --custom_decontamination_index contaminants_index
  ```
- Faster (less sensitive) alignment settings: `--bowtie_k 1 --bowtie_score_min "L,-0.6,-0.6"`

---

## Cloud Execution Issues

### AWS: Access denied

**Error:**
```
Access Denied (Service: Amazon S3)
```

**Solution:**
1. Check AWS credentials:
   ```bash
   aws sts get-caller-identity
   ```

2. Verify S3 bucket permissions
3. Check IAM role attached to Batch compute environment

### AWS: Batch job failed

**Error:**
```
Essential container in task exited
```

**Solution:**
1. Check CloudWatch logs for the job
2. Verify compute environment has sufficient resources
3. Check container can access S3 paths

### GCP: Quota exceeded

**Error:**
```
Quota exceeded for resource
```

**Solution:**
1. Request quota increase in GCP Console
2. Use smaller instance types
3. Reduce parallelism with `-queue-size`

### Azure: Authentication failed

**Error:**
```
Azure Batch authentication failed
```

**Solution:**
1. Verify environment variables:
   ```bash
   echo $AZURE_BATCH_ACCOUNT_NAME
   echo $AZURE_STORAGE_ACCOUNT_NAME
   ```

2. Check account keys are correct
3. Verify account is in correct region

---

## Output Issues

### Missing output files

**Symptom:** Expected output files not present

**Solution:**
1. Check if process completed:
   ```bash
   cat .nextflow.log | grep -i error
   ```

2. Check work directory for intermediate files

### Corrupted output files

**Symptom:** Output files have zero size or are truncated

**Solution:**
1. Check disk space during execution
2. Verify process exit status in trace file
3. Re-run with `-resume`

---

## Debug Mode

### Enable verbose logging

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --output ./results \
    -with-trace \
    -with-report \
    -with-timeline \
    -with-dag \
    -dump-channels
```

### Inspect failed process

```bash
# Find work directory of failed task
cat .nextflow.log | grep -A5 "Error executing process"

# Inspect work directory
ls -la work/xx/xxxxxxxx/
cat work/xx/xxxxxxxx/.command.log
cat work/xx/xxxxxxxx/.command.err
```

### Test with stub mode

```bash
# Dry run without executing actual commands (test profile supplies the input/output params)
nextflow run main.nf -profile test,docker -stub
```

### Preview the workflow graph

```bash
# Resolve parameters and wiring without executing any process
nextflow run main.nf -profile docker --input samplesheet.csv --output ./results -preview
```

---

## Quick Fixes Checklist

- [ ] Using the latest pipeline version
- [ ] Samplesheet exists and is formatted correctly (`sample,r1,r2,s` header)
- [ ] FASTA/FASTQ files exist and paths are absolute
- [ ] Sufficient disk space and memory allocated
- [ ] YAML syntax is correct (comma-separated strings, no spaces after commas)
- [ ] Work directory exists (for `-resume`)
- [ ] Container runtime works (Docker/Singularity)

---

## Getting Help

### Collect debug information

Before asking for help, gather:

1. **Nextflow version:** `nextflow -version`
2. **Error message:** Copy full error output
3. **Log file:** `.nextflow.log`
4. **Trace file:** `pipeline_info/execution_trace_*.txt`
5. **Command used:** Full nextflow command

### Resources

- **GitHub Issues:** https://github.com/gene2dis/BugBuster/issues
- **Nextflow Documentation:** https://nextflow.io/docs/latest/
- **Nextflow Slack:** https://www.nextflow.io/slack-invite.html
- **Contact:** ffuentessantander@gmail.com

### Reporting bugs

When reporting bugs, include:

1. Minimal reproducible example
2. Expected vs actual behavior
3. Environment details (OS, container runtime, cloud provider)
4. Relevant log snippets
