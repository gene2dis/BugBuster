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
ERROR: Nextflow version 24.04.0 or later is required
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
- If you are following instructions written for an older release, the parameter may have been removed or renamed (e.g. `--kraken_db_used`, `--sourmash_db_name`, `--validationShowHiddenParams`, and all `--mmseqs_*` parameters no longer exist; `--bbmap_lenght` is now `--bbmap_length`)

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
