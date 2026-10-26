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

### `--<option> false` is ignored, or a launch check fires for an option you turned off

**Symptoms (older pipeline versions on Nextflow 26.04 or newer):** a branch runs although
you passed `--<option> false`, or the launch fails with a dependency check such as
`--include_binning requires an assembly` for an option you set to `false`; a numeric option
such as `--min_read_sample 1000` fails with `Cannot compare java.lang.Integer ... with
java.lang.String`.

**Cause:** since Nextflow 26.04 every command-line parameter reaches the pipeline as text
(`false` is the string `"false"`), and the parameter schema check does not convert it. Nextflow
25.10 still converted these values itself.

**Solution:** fixed in the current release, which reads every on/off option and numeric
comparison in a way that works on both Nextflow versions. On an older release, omit the flag
instead of passing `false` (all on/off options except `--quality_control`, `--microbecensus`,
`--functional_cazy`, `--rgi_include_wildcard` and `--azure_delete_pools` are off by default),
or set the values in a `-params-file` YAML, where `false` and numbers keep their types.

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

   Co-assembly MEGAHIT already requests 6 GB per GB of gzipped pooled reads
   (at least 72 GB), but every request is capped by `--max_memory` (default
   128 GB): for a large co-assembly, raise it to the node's RAM first.
   `--megahit_min_count 3` lowers MEGAHIT's graph-build memory (17 % in the
   assembly benchmark, for 5 % less assembled sequence), and
   `--assembly_read_reduction normalize` shrinks the assembler input.

2. For specific processes, put the setting in a config file and pass it with
   `-c custom.config` (requests above `--max_memory` are still capped by it):
   ```groovy
   process {
       withName: 'MEGAHIT' {
           memory = '128.GB'
       }
   }
   ```

3. Use a larger instance type (cloud) or node (HPC)

A process that runs out of memory is retried twice with more memory each time
(exit statuses 130-145 and 104). MEGAHIT reports running out of memory as 137
too (`bin/megahit_run.sh`), and its graph build is held to 90% of the task
memory. MEGAHIT takes the minimum a step needs even above that budget, so a
deep co-assembly can still need a larger request: its log
(`<work dir>/coassembly_assembly/log`) gives each step's "Minimum memory
required" in bytes.

The heaviest functional-annotation steps are `EGGNOG_MAPPER_SEARCH` (~80 GB peak
on a real co-assembly; it requests 96 GB on the first attempt), `WOLTKA_ALIGN` (≥ 68 GB,
see below) and `HUMANN` (~23.5 GB per sample). Observed requirements per step are
listed under [Functional annotation branches](manual.md#functional-annotation-branches)
in the manual.

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

3. Run fewer tasks at once (lower `--max_cpus` with the local executor) or split the
   samples into batches; see [`DISK_OPTIMIZATION.md`](DISK_OPTIMIZATION.md#running-with-limited-disk-space).

Note that `cleanup = true` (the `low_disk` profile) does not help with this error: it
deletes the work directory only after the run succeeds, so peak usage during the run is
unchanged.

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

With the local executor, a task killed at its time limit can instead fail with
`process hasn't exited` (`java.lang.IllegalThreadStateException` in `.nextflow.log`), and its
`.exitcode` is 143 (written after the kill). Nextflow then has no exit status for the task, so
the error strategy stops the run instead of retrying it with more time. Compare the task's
start and end times with its `time` setting; if they match, raise `time` for that process as
above and `-resume`.

---

### "Process requirement exceeds available CPUs" (or memory), or raised limits seem ignored

**Error:**
```
Process requirement exceeds available CPUs -- req: 8; avail: 2
```

**Cause:** With the local executor, `--max_cpus`/`--max_memory` set two things: the cap on
every task, and the run's machine-wide CPU and memory pool. The pool is fixed before a `-c`
config file is read. A `-c` file that raises `max_*` in `params {}` lifts the per-task caps
but leaves the pool at the profile's values, so a task bigger than the pool cannot be
scheduled, and the run stops at the first such task. After `-profile test`, the pool is
2 CPUs / 6 GB.

**Solution:** either pass the limits on the command line or in a `-params-file` (these set
both), or add an `executor` block to the `-c` file:
```bash
nextflow run main.nf -profile test,docker --max_cpus 32 --max_memory 256.GB --max_time 48.h ...
```
```groovy
// my_config.config, passed with -c
params   { max_cpus = 32; max_memory = '256.GB'; max_time = '48.h' }
executor { cpus = 32; memory = '256.GB' }
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

### eggNOG-mapper image cannot be pulled under docker/podman

**Error:** `EGGNOG_MAPPER_SEARCH` or `EGGNOG_MAPPER_ANNOTATE` fails before
running, with `denied`, `unauthorized` or `manifest unknown` for
`ghcr.io/gene2dis/bugbuster-eggnog-mapper`.

**Cause:** eggNOG-mapper 3.0.0-beta6 has no official docker image, so under
docker/podman (and the docker-based cloud profiles) the contig branch uses an
image of the same version built by this pipeline and hosted on GitHub's
container registry. The pull fails if the machine cannot reach `ghcr.io`, or if
a registry login for another account is cached and rejected.

**Solution:** check that `docker pull ghcr.io/gene2dis/bugbuster-eggnog-mapper@<digest>`
(digest from `modules/local/eggnog_mapper_search/main.nf`) works on the machine
running the tasks; the image is public, so no login is needed (`docker logout ghcr.io`
clears a stale one). Offline hosts can load the image from a saved archive
(`docker save` / `docker load`). Alternatively run with `-profile apptainer` or
`-profile singularity`, which use the official eggNOG-mapper Apptainer image.

> Related: the official eggNOG-mapper Apptainer image ships without `procps`,
> which would normally abort Nextflow tasks with `Command 'ps' required by
> nextflow to collect task metrics cannot be found`. The pipeline works around
> this with a delegating `ps` shim in `bin/` — do not remove `bin/ps` while the
> singularity/apptainer path uses that image (the docker image includes `ps`).

### EGGNOG_MAPPER_ANNOTATE fails with "the installed eggNOG-mapper _parse_ogs_string is not the v3.0.0-beta6 method"

**Cause:** the annotate step runs eggNOG-mapper through `bin/emapper_ogs_fix.py`,
which corrects a bug in eggNOG-mapper 3.0.0-beta6 (upstream issue #620): OG
names whose family part contains `|` (about 2 % of eggNOG 7 OGs, e.g.
`ABC_tran|TL31Y9@131567|A-1`) were cut short, so those genes lost their COG
category or took it from another OG (about 1-4 % of annotated genes; KO, EC,
Pfam and CAZy were not affected). The wrapper replaces only that parsing
function, and only after checking that the installed code is exactly the
beta6 version it was written for. This error means the eggNOG-mapper image
was changed to a different version.

**Solution:** this only happens if the eggNOG-mapper container was changed
(e.g. a custom `process.container` override, or a pipeline update to a newer
eggNOG-mapper). Use the pinned image, or — when updating eggNOG-mapper —
check whether issue #620 is fixed in the new version and update or remove
`bin/emapper_ogs_fix.py` accordingly.

### RGI_BWT fails: "KMA failed", "bowtie2 failed" or "not the 6.0.8 method"

`RGI_BWT` runs `rgi bwt` through `bin/rgi_bwt_kma_fix.py`, which keeps only the reads the aligner placed on
a reference, for KMA grouped by read name (RGI keeps every read: tens of GB of disk per deep sample and the
whole read table in memory), and checks that the aligner succeeded. Its messages:

- **"KMA failed (exit status N)":** the alignment failed (status 2: a missing or incomplete KMA index in
  `localDB/`; 137: killed, usually out of memory, retried with more memory). RGI itself ignored this and
  reported the sample as having no hits. KMA 1.4.9 (RGI 6.0.3's image) could exit 95 after a successful
  run, so 95 is still accepted when KMA has written its result table.
- **"bowtie2 failed" / "bwa failed (exit status N)"** (`--rgi_aligner bowtie2|bwa`): the alignment
  failed, usually a missing or incomplete index in the task's `localDB/`.
- **"filtering the KMA SAM failed"** (or the bowtie2/bwa SAM): the filter, or for KMA the sort that groups
  the reads by name, failed; awk exit status 137 means it was killed, usually out of memory, and the task is
  retried with more memory.
- **"the installed RGI BWT.<method> is not the 6.0.8 method this patch targets":** the RGI container was
  changed. Use the pinned image, or, when updating RGI, check whether the new version still writes every read
  to the SAM, then update or remove `bin/rgi_bwt_kma_fix.py`.

### RGI_KMER fails: "worker was killed", "did not finish within", or "not the 6.0.8 method"

`RGI_KMER` runs `rgi kmer_query` through `bin/rgi_kmer_kma_fix.py`. That wrapper
fixes an RGI bug that made the task hang forever: the run waited for a worker
that had died. (It also fixed RGI 6.0.3's reading of KMA BAMs, which RGI 6.0.8
does correctly.) It also stops a run that takes too long. Its messages:

- **"an rgi kmer_query worker was killed by signal 9 (out of memory?)" (exit 137):**
  a worker ran out of memory. Each worker is a copy of the process holding the
  ~17 GB k-mer database, so memory grows with the number of workers. Nextflow
  retries the task with more memory (36 GB, then 72 GB, then 108 GB) automatically.
  If every attempt fails, raise the memory, or keep the worker count low (the
  default is 2 CPUs = 2 workers):

  ```groovy
  // custom.config, used with -c custom.config
  process {
      withName: 'RGI_KMER' {
          memory = 96.GB
      }
  }
  ```

- **"rgi did not finish within 7200 s" (exit 124):** the run hit the 2 h limit
  (not retried). On the test data `kmer_query` takes about 3.5 minutes; a very
  large sample can take longer. Raise the limit in a config file:

  ```groovy
  process {
      withName: 'RGI_KMER' {
          ext.time_limit = 8.h
      }
  }
  ```

- **"the installed RGI CARDkmers.<method> is not the 6.0.8 method this patch
  targets":** the RGI container was changed (a `process.container` override, or a
  pipeline update to a newer RGI). Use the pinned image, or, when updating RGI,
  check whether the new version still needs the fixes, then update or remove
  `bin/rgi_kmer_kma_fix.py`.

### KMARVEL fails: "protein BLAST step failed", "did not finish", or "has no CARD_protein_all.fasta"

`KMARVEL` runs K-MARVEL and then `bin/kmarvel_check_output.sh` on its log,
because K-MARVEL 1.0.0 exits 0 even when its blastp step fails.

- **"K-MARVEL found N ARGs but its protein BLAST step failed":** blastp could
  not run on the assembled ARGs, so the table would be unfiltered and lack the
  `blast_*` columns. The task fails instead of publishing it. See
  `{sample}_kmarvel.log` in the task work directory; out-of-memory or a broken
  BLAST index of `CARD_protein_all.fasta` are the usual causes.
- **"no 'Found N Antimicrobial Resistance Genes' line ... K-MARVEL did not
  finish":** K-MARVEL stopped before reporting its ARG count (a Rust panic, for
  example an unexpected character in the reads), or a different K-MARVEL
  version changed its log. Check `.command.log` in the task work directory.
- **"... has no CARD_protein_all.fasta (expected a FORMAT_KMARVEL_DB /
  --custom_kmarvel_db layout)":** `--custom_kmarvel_db` must point at the
  directory holding `CARD_protein_all.fasta`, `metadata.tsv` and `snps.txt`
  (for example `databases/kmarvel/kmarvel_db` from an earlier run), not at its
  parent.
- **No ARGs for any sample:** K-MARVEL ignores reads shorter than
  `--kmarvel_min_read_length` (100 bp), so heavily trimmed short reads can leave
  little to assemble. Lowering `--kmarvel_min_coverage` / `--kmarvel_min_identity`
  (90 % each) admits partial or more divergent genes.

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

Two related guards protect the term parsing (design: never misparse silently):
`PFAMs value '<X>' is not '<pfam_name>_<start>_<end>'` and
`COG_category value '<X>' is neither category letters nor a COG id` mean the
eggNOG-mapper output format changed — re-verify before extending the parser.
`COG id '<X>' is not in the --cog-def table` means the COG definitions table
(`--cog_db` / `--custom_cog_db`) is older than the eggNOG database's COG ids:
use a COG release that covers them (the default `cog-24` covers all of eggNOG 7).

The Woltka read backend has the matching guard in `WOLTKA_CLASSIFY`:
`COG id '<X>' (<KO> in function/kegg/ko-to-cog.map) is not in the COG definitions
table ... and is not a known WoLr2 defect` means the WoL database and the COG
table do not match — the default `cog-24` covers every WoLr2 COG id except nine
known defects of that release, which are skipped with a `WARNING: ... known WoLr2
defects` line in the task log (expected, not an error). A different WoL release or
COG table must be re-checked together before extending the skip-list
(`WOLR2_UNMAPPABLE_COGS` in `bin/woltka_function_profile.py`). `Pfam accession
<X> ... has no entry in function/pfam/pfam_name.txt` and `Pfam name <X> is shared
by ...` mean the WoL Pfam maps are incomplete or from a mismatched release.

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
--mag_level_functional requires at least two --binners: Binette refinement and its completeness/contamination quality filter only run with >= 2 binners ...
```

**Cause:** deliberate validation, not a bug. Bakta annotates the refined bins,
so binning must be enabled — and the pipeline's only bin quality filter is
Binette refinement (`--bin_min_completeness` / `--bin_max_contamination`,
defaults 50/10), which only runs when at least
two binners are selected. With a single binner the bins channel carries raw,
unfiltered binner output, and the MAG branch refuses to annotate unfiltered
bins.

**Solution:** enable binning with two or more binners:
```bash
nextflow run main.nf ... --include_binning --binners semibin,metabat2 --mag_level_functional
```
Note that Bakta annotates every Binette-refined bin; it does not further
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
    quay.io/biocontainers/bakta:1.12.2--pyhdfd78af_0 \
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
off entirely. The read-level branch (`--read_level_functional`) shares the same
MicrobeCensus output: there the affected sample has an empty
`abundance_cpge` in `read_function_abundance.tsv` and
`cpge_status = unavailable` in `read_sample_summary.tsv`, while its native read
counts are complete.

**Note:** MicrobeCensus is invoked through `bin/run_microbe_census_py3fix.py`.
Both published biocontainers of the unmaintained upstream tool are broken
(the Python-3 image lacks an unreleased upstream str/bytes fix — upstream
issue #32 — and the Python-2 image cannot load its bundled RAPsearch2
binary); the shim patches the one broken function in the pinned Python-3
image and delegates everything else to the stock tool. Do not remove the shim
while the pin is `microbecensus:1.1.1--pyhca03a8a_2`.

### WOLTKA_ALIGN runs out of memory / is killed

**Symptom:** `WOLTKA_ALIGN` fails with exit status 137/140 or a Bowtie2
"Out of memory" message while loading the index.

**Cause:** the Web of Life (WoLr2) Bowtie2 index is 93.6 GB on disk and needs
**at least 68 GB RAM** to load. The task requests up to 256 GB (labels
`process_high` + `process_high_memory`) but is capped by `--max_memory`.

**Solution:** run on a node with ≥ 68 GB available and set `--max_memory`
accordingly (e.g. `--max_memory 96.GB`). There is no smaller WoLr2 index for
this backend.

### AGGREGATE_READ_FUNCTIONS / WOLTKA_CLASSIFY fail with a layout or database message

**Errors:**
```
woltka version '<X>' is not among the versions this parser was verified against ...
... unexpected Woltka profile header ...
<N> profile ORF(s) have no entry in proteins/length.map.xz ...
Woltka database file missing: .../function/...
```

**Cause:** deliberate guards (design: silent misparsing is the failure mode
prevented). The first two mean the Woltka version or its output layout is not
the one the read-branch scripts were verified against (pinned: woltka 0.1.7) —
re-verify on real output before extending `KNOWN_WOLTKA_VERSIONS` in
`bin/aggregate_read_functions.py`, plus `tests/bin/test_aggregate_read_functions.sh`.
The last two mean `--custom_woltka_db` is incomplete or mixes releases: the
Bowtie2 index and `proteins/`/`function/` files must all come from the same
WoLr2 release, in the FTP layout.

**Solution:** for database errors, complete the mirror with the file list in
`docs/parameters.md` (`--custom_woltka_db`) or the download recipe in
`docs/manual.md`, and check the md5 files (`xz -dc <file>.xz | md5sum` against
`<file>.md5` — WoLr2 checksums cover the uncompressed content).

### SUPERFOCUS fails: not a SUPER-FOCUS database root / no DB_90 files for the aligner

**Errors:**
```
ERROR: <dir> looks like the db/ folder itself; point --custom_superfocus_db at its PARENT (the directory containing db/)
ERROR: <dir>/db/database_PKs.txt not found - not a SUPER-FOCUS database root
ERROR: <dir>/db/static/mmseqs2/90_clusters.db not found - the database has no mmseqs2 DB_90 files (--superfocus_aligner mmseqs2)
```

**Cause:** `--custom_superfocus_db` must be the database **root**, the directory
that contains `db/` (`db/database_PKs.txt` plus `db/static/<aligner>/`), not
`db/` itself. The database formats are aligner-specific, so the root must also
hold the DB_90 files of the selected `--superfocus_aligner`
(`db/static/diamond/90_clusters.db.dmnd` or `db/static/mmseqs2/90_clusters.db*`).

**Solution:** pass the parent of `db/`; if the aligner's files are missing,
either switch `--superfocus_aligner` to the one the database was built for or
add that aligner's DB_90 archive (download recipe in `docs/manual.md`, Manual
Database Download).

### AGGREGATE_READ_FUNCTIONS / SUPERFOCUS fail with a layout or version message

**Errors:**
```
superfocus version '<X>' is not among the versions this parser was verified against ...
... unexpected header ... the SUPER-FOCUS output layout changed or the query was not the single file '<sample>.fastq'
... total count <N> exceeds the <M> input reads - was the table produced with -n 0, or for a different query?
```

**Cause:** deliberate guards, as for Woltka. The first two mean the SUPER-FOCUS
version or its output layout is not the one the scripts were verified against
(pinned: SUPER-FOCUS 1.8) — re-verify on real output before extending
`KNOWN_SUPERFOCUS_VERSIONS` in `bin/aggregate_read_functions.py` (and the layout
checks in `bin/superfocus_function_profile.py`, plus
`tests/bin/test_superfocus_scripts.sh`). The last means a read was counted more
than once: the per-level composition relies on SUPER-FOCUS dividing each read
1/k (`-n 1`, fixed in the module); an `-n 0` added to the `SUPERFOCUS`
`ext.args` breaks that.

**Solution:** keep the pinned container and remove `-n` from any custom
`ext.args` for `SUPERFOCUS` (only the identity / alignment-length / e-value /
fast-mode thresholds belong there).

### HUMANN fails: database root, MetaPhlAn database version, or sample name

**Errors:**
```
ERROR: <dir>/metaphlan/ not found - not a HUMAnN database root (chocophlan/, uniref/, utility_mapping/, metaphlan/; see --custom_humann_db)
ERROR: <dir>/metaphlan/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl not found - HUMAnN 4.0.0a2 accepts only the MetaPhlAn database mpa_vOct22_CHOCOPhlAnSGB_202403
ERROR: MetaPhlAn DB version check failed. Expected one of: ['vOct22_CHOCOPhlAnSGB_202403'] Detected tag-like strings: [...]
Invalid sample name '<name>' for --read_level_functional humann: sample names must not contain 's__' or 't__' ...
```

**Cause:** `--custom_humann_db` must be a database **root** with `chocophlan/`,
`uniref/`, `utility_mapping/` and `metaphlan/`. HUMAnN 4.0.0a2 (the pinned
alpha) accepts exactly one MetaPhlAn database, `mpa_vOct22_CHOCOPhlAnSGB_202403`,
and checks its tag at run time — a MetaPhlAn 4 database you already have (e.g.
`vJun23`, or the server's newer `mpa_latest`) is refused. The pipeline exposes
`<root>/metaphlan/` to that check through `METAPHLAN_DB_DIR`; the files must sit
directly in that folder. The sample-name rule works around a 4.0.0a2 defect: a
MetaPhlAn profile line holding both `s__` and `t__` — including the header line
that carries the sample name — is misread as a taxon row and crashes HUMAnN with
`IndexError: list index out of range`.

**Solution:** let the pipeline download the database set (`--humann_db`), or
assemble the root as described for `--custom_humann_db` in `docs/parameters.md`
(the MetaPhlAn files come from
`https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/`:
`mpa_vOct22_CHOCOPhlAnSGB_202403.tar` and
`bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar`, unpacked flat into
`metaphlan/`). Rename samples whose names contain `s__` or `t__`.

### AGGREGATE_READ_FUNCTIONS / HUMANN fail with a layout or version message

**Errors:**
```
... unexpected gene-family header '...' - not the HUMAnN layout this parser was verified against ...
... HUMAnN version '<X>' is not among the versions this parser was verified against ...
... unexpected gene-family feature / unexpected pathway feature ...
humann version '<X>' is not among the versions this parser was verified against ...
```

**Cause:** deliberate guards (design doc Section 4.6.1: HUMAnN 4 is an alpha
whose output layout may change between builds; silent misparsing is the failure
being prevented). The pinned build is HUMAnN 4.0.0a2, run with
`--count-normalization RPKs`; a header that names other units (e.g. `Adjusted
CPMs`) means `--count-normalization` was overridden in the `HUMANN` `ext.args`.

**Solution:** keep the pinned container and do not set `--count-normalization`
in `ext.args`. Moving to another HUMAnN build requires re-verifying its output
on real data and then extending `KNOWN_HUMANN_VERSIONS` in both
`bin/humann_function_profile.py` and `bin/aggregate_read_functions.py` (plus
`tests/bin/test_humann_scripts.sh`). Note also that HUMAnN 4.0.0a2's own
`humann_renorm_table` / `humann_regroup_table` mishandle the `READS_UNMAPPED`
row of its gene-family table (fixed upstream after 4.0.0a2): if you post-process
the raw `_2_genefamilies.tsv` yourself, drop that row first.

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
       --custom_sylph_db /path/to/gtdb-r232-c200-dbv2.syl2db \
       --custom_sylph_taxonomy /path/to/gtdb_r232_metadata.tsv.gz \
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

### Sylph: Classified is 1.0, or TAXONOMY_PHYLOSEQ fails with "No sylph genomes in any sample"

**Symptom (a):** with `--taxonomic_profiler sylph`, the `Classified` column of
`02_taxonomy/tables/Reads_report_taxonomy.csv` is 1.0 (or close) for samples that surely hold organisms
missing from the database, and the estimated read pairs in
`sylph_<db>_read_counts.tsv` add up to every read pair.

**Cause:** the classified fraction is the sum of sylph's `Sequence_abundance`
under `-u`, which rests on sylph's estimate of the read accuracy. Reads from
genomes absent from the database can lower that estimate, so sylph
overestimates coverage and caps the explained fraction at 100% (seen on
simulated data; sylph itself warns that `-u` has caveats for complex
communities). The relative abundance table is not affected.

**Solution:** set the read identity explicitly with `--sylph_read_seq_id`
(Illumina: about 99.5) and rerun with `-resume`, which re-runs only the sylph
profiling and the steps after it.

**Symptom (b):** `TAXONOMY_PHYLOSEQ` fails with *No sylph genomes in any
sample*. No genome passed sylph's detection threshold in any sample (each
skipped sample is warned about first): very shallow samples, or a database that
does not cover the community. Check the per-sample `.tsv` profiles in
`02_taxonomy/sylph/` and the `Classified` column, and use a database that fits
the samples (`--sylph_db`, or `--custom_sylph_db` together with its
`--custom_sylph_taxonomy`), or rerun with `--taxonomic_profiler none`.

---

## Decontamination Issues

The default aligner is bwa-mem3 (`BWAMEM3_*` tasks); `--read_aligner bowtie2`
runs the `BOWTIE2_*` tasks instead.

### A BWAMEM3_* task fails with "Illegal instruction" (exit status 132)

The bioconda bwa-mem3 build needs an x86-64 CPU with AVX2 (Haswell 2013 or
newer). Check a compute node with `grep -c avx2 /proc/cpuinfo` (0 = no AVX2)
and run with `--read_aligner bowtie2` on nodes without it.

### Process BWAMEM3_BUILD_COMBINED failed

**Exit status 3 / "needs more memory":** the human index build needs ~75 GB
RAM; the pipeline warns at launch when `--max_memory` is lower. Either:

- raise `--max_memory` (the build runs once; later runs reuse the index with
  `--custom_decontamination_index <databases_dir>/bwamem3_index/contaminants_index`),
- pass an index built elsewhere (`bwa-mem3 index -p contaminants_index/contaminants contaminants.fasta`,
  or a bwa-mem2 index), or
- use `--read_aligner bowtie2` (~6 GB to build, about an hour).

For other failures, check the FASTA files as for BOWTIE2_BUILD_COMBINED below.

### BWAMEM3_DECONTAMINATE fails with "no bwa-mem3/bwa-mem2 index", or a Bowtie2 index is rejected at launch

`--custom_decontamination_index` points to a Bowtie2 index (`*.bt2`), for
example one built by an earlier release, where bowtie2 was the default. Add
`--read_aligner bowtie2` to keep using it, or point the option at a bwa-mem3
index (`*.bwt.2bit.64`).

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
   # Should list contaminants.1.bt2 ... contaminants.rev.2.bt2 (.bt2l for a large index)
   # (the default bwa-mem3 index is in bwamem3_index/contaminants_index/:
   #  contaminants.amb .ann .bwt.2bit.64 .pac)
   ```
2. Verify read files decompress cleanly: `zcat sample_R1.fastq.gz | head -n 4`
3. Check disk space: `df -h .`

### Multiple host genomes

`custom_host_fasta` takes a **comma-separated string**, not a YAML list; spaces
around the commas are ignored, and the paths themselves may not contain commas
or spaces. Earlier versions read the whole string as one path, so a list failed
with a missing-file error.

```yaml
# Correct - comma-separated string
custom_host_fasta: "/path/file1.fasta,/path/file2.fasta,/path/file3.fasta"
custom_host_fasta: "/path/file1.fasta, /path/file2.fasta"

# Incorrect - YAML list (not supported)
custom_host_fasta:
  - "/path/file1.fasta"
  - "/path/file2.fasta"
```

All the genomes go into one index; the bwa-mem3 build requests ~24 GB RAM per
Gbp of combined genome, so raise `--max_memory` accordingly.

### Slow index building or decontamination

- Pre-build the index once and reuse it:
  ```bash
  cat phix.fasta host.fasta > contaminants.fasta
  mkdir -p contaminants_index
  bwa-mem3 index -t 16 -p contaminants_index/contaminants contaminants.fasta

  nextflow run main.nf ... --custom_decontamination_index contaminants_index
  ```
- Under `--read_aligner bowtie2`, the human index build takes about an hour
  (`bowtie2-build --threads 16 contaminants.fasta contaminants_index/contaminants`;
  pass the index with `--read_aligner bowtie2`). bowtie2 runs with `--local`, so
  `--bowtie_score_min` takes the local form `G,<a>,<b>` (default `G,15,6`)

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

### No bins: empty bin tables, `SKIPPED.txt` instead of bins

**Symptom:** the run succeeds, but `04_binning/` has no bins. A bin folder holds
only `SKIPPED.txt`. The bin tables (`Mag_quality_summary.csv`,
`MAGs_tax_summary.csv`, `Bin_summary.csv`) have a header and no rows, the plots
say "No bins were assessed" / "No MAGs were classified", and the log shows
warnings such as `all CheckM2 quality reports are empty` or
`No refined bins for coassembly: skipping bin depth`.

**Cause:** binning found nothing usable on this data. This is common for
host-dominated or shallow samples, where little microbial DNA is left after
host removal:

- `SemiBin found no bins`: SemiBin2 had no bin before reclustering.
- `All contigs < <N>bp - SemiBin skipped`: no contig is long enough for SemiBin2.
- `No bin reached <completeness>% completeness and <contamination>% contamination`:
  Binette ran, but no refined bin passed `--bin_min_completeness` /
  `--bin_max_contamination` (default 50 / 10).
- `No bins available for refinement`: no binner produced a bin for Binette.

The run finishes on purpose instead of failing. A report step still fails when
an input file is missing, which points to a real pipeline problem.

**What to try:** check the clean-read counts in `Reads_report.csv` and the
assembly size first. Co-assembling samples (`--assembly_mode coassembly`) or
sequencing deeper gives the binners more coverage. Lowering the refinement
thresholds keeps lower-quality bins, at the cost of more incomplete or
contaminated MAGs.

### Corrupted output files

**Symptom:** Output files have zero size or are truncated

**Solution:**
1. Check disk space during execution
2. Verify process exit status in trace file
3. Re-run with `-resume`

### `-resume` re-runs a summary or batch step that already finished

**Symptom (older pipeline versions):** a resumed run that changed nothing still re-runs a step
that combines many samples. The steps affected were the co-assembly (`MEGAHIT` / its read
alignment), `METAWRAP`, `CHECKM2_BATCH`, `GTDB_TK_BATCH`, or one of the report steps (reads,
taxonomy, ARG, RGI and bin reports).

**Cause:** these steps received their inputs in the order the upstream tasks happened to
finish, which can differ between runs and changes the step's cache key.

**Solution:** fixed in the current release; inputs are now passed in a fixed order (by sample).
After upgrading, an existing run may re-run each of these steps once on its next `-resume`,
then stays fully cached. Note that MEGAHIT itself is not byte-reproducible with several
threads (a re-run can differ in a handful of contigs), so a co-assembly that re-runs is not
expected to be identical to the previous one.

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
- [ ] YAML syntax is correct (lists such as `custom_host_fasta` and `binners` are comma-separated strings)
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
