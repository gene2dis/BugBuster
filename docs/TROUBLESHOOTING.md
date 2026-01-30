# Troubleshooting Guide

## Common Errors and Solutions

### Error: "Cannot emit a multi-channel output: decontamination_index"

**Error Message:**
```
ERROR ~ Cannot emit a multi-channel output: decontamination_index
java.lang.IllegalArgumentException: Cannot emit a multi-channel output: decontamination_index
```

**Cause:** This error occurred in earlier versions where the `decontamination_index` channel was being assigned in multiple conditional branches.

**Solution:** This has been fixed in the latest version. Update your pipeline:
```bash
cd /path/to/BugBuster
git pull
```

If you're still experiencing this issue, ensure you're using the latest version of the PREPARE_DATABASES subworkflow.

---

### Error: "Cannot find samplesheet file"

**Error Message:**
```
Cannot find samplesheet file: /path/to/samplesheet.csv
```

**Cause:** The samplesheet file doesn't exist at the specified path.

**Solution:**
1. Verify the file exists:
   ```bash
   ls -lh /path/to/samplesheet.csv
   ```

2. Use absolute path:
   ```bash
   nextflow run main.nf --input /absolute/path/to/samplesheet.csv
   ```

3. Check samplesheet format:
   ```csv
   sample,fastq_1,fastq_2
   sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
   sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
   ```

---

### Error: "Process BOWTIE2_BUILD_COMBINED failed"

**Possible Causes:**
1. FASTA files not found
2. Insufficient memory
3. Corrupted FASTA files

**Solutions:**

**1. Verify FASTA files exist:**
```bash
# Check phiX FASTA
ls -lh /path/to/phix.fasta

# Check host FASTA
ls -lh /path/to/host.fasta
```

**2. Increase memory allocation:**
```bash
nextflow run main.nf \
  --input samples.csv \
  --output results \
  --max_memory 64.GB
```

**3. Validate FASTA files:**
```bash
# Check if file is valid
zcat phix.fasta.gz | head -n 10

# Test decompression
gunzip -t phix.fasta.gz
```

---

### Error: "Process BOWTIE2_DECONTAMINATE failed"

**Possible Causes:**
1. Index not built correctly
2. Corrupted read files
3. Insufficient disk space

**Solutions:**

**1. Check index files:**
```bash
ls -lh databases/bowtie_index/contaminants_index/
# Should see: contaminants.1.bt2, contaminants.2.bt2, etc.
```

**2. Verify read files:**
```bash
zcat sample_R1.fastq.gz | head -n 4
```

**3. Check disk space:**
```bash
df -h .
```

---

### Error: High Memory Usage

**Symptom:** Process killed due to memory limits

**Solution:**

**For index building:**
```bash
nextflow run main.nf --max_memory 128.GB
```

**For specific process:**
Edit `nextflow.config`:
```groovy
process {
    withName: 'BOWTIE2_BUILD_COMBINED' {
        memory = '64.GB'
    }
    withName: 'BOWTIE2_DECONTAMINATE' {
        memory = '32.GB'
    }
}
```

---

### Error: "Channel has been used as output by more than one operator"

**Cause:** Channel reuse issue in workflow

**Solution:** Use `.first()` or create separate channels:
```groovy
// Instead of reusing channel
ch_index = some_process.out.index

// Use first() for single-value channels
ch_index = some_process.out.index.first()
```

---

### Warning: "The use of Channel to access channel factories is deprecated"

**Warning Message:**
```
The use of `Channel` to access channel factories is deprecated -- use `channel` instead
```

**Cause:** Using uppercase `Channel` instead of lowercase `channel`

**Solution:** This is a cosmetic warning and doesn't affect functionality. It will be addressed in future updates.

---

### Error: Multiple Genomes Not Being Used

**Symptom:** Only one genome is used for decontamination

**Solution:**

**Check YAML format:**
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

**Verify index building:**
```bash
# Check work directory for concatenated FASTA
find work -name "combined_contaminants.fasta" -exec wc -l {} \;
```

---

### Error: Resume Not Working

**Symptom:** Pipeline restarts from beginning with `-resume`

**Possible Causes:**
1. Work directory deleted
2. Configuration changed
3. Input files modified

**Solutions:**

**1. Verify work directory exists:**
```bash
ls -lh work/
```

**2. Use same configuration:**
```bash
# Use same parameters and profile
nextflow run main.nf -params-file params.yaml -profile docker -resume
```

**3. Check cache:**
```bash
# View cached tasks
nextflow log last
```

---

### Performance Issues

#### Slow Index Building

**Solution 1: Pre-build index**
```bash
cat phix.fasta host.fasta > contaminants.fasta
bowtie2-build --threads 16 contaminants.fasta contaminants_index/contaminants

nextflow run main.nf --custom_decontamination_index contaminants_index
```

**Solution 2: Increase CPUs**
```bash
nextflow run main.nf --max_cpus 32
```

#### Slow Decontamination

**Solution: Adjust Bowtie2 parameters**
```yaml
# In params.yaml - faster but less sensitive
bowtie_k: 1
bowtie_score_min: "L,-0.6,-0.6"
```

#### High Disk Usage

**Solution: Enable cleanup**
```bash
nextflow run main.nf -profile low_disk
```

Or in `nextflow.config`:
```groovy
params {
    enable_work_cleanup = true
    store_clean_reads = true
}
```

---

## Debugging Tips

### View Detailed Logs

```bash
# View last run log
nextflow log last

# View specific run
nextflow log <run_name>

# View with details
nextflow log last -f name,status,duration,realtime
```

### Check Process Status

```bash
# View execution report
open results/pipeline_info/execution_report_*.html

# View timeline
open results/pipeline_info/execution_timeline_*.html

# View resource usage
cat results/pipeline_info/execution_trace_*.txt
```

### Inspect Work Directory

```bash
# Find specific process work directory
find work -name "*.command.sh" | grep BOWTIE2_BUILD_COMBINED

# View process script
cat work/xx/xxxxx/.command.sh

# View process output
cat work/xx/xxxxx/.command.out

# View process error
cat work/xx/xxxxx/.command.err
```

### Test Configuration

```bash
# Dry run
nextflow run main.nf -params-file params.yaml -preview

# Validate YAML
python3 -c "import yaml; print(yaml.safe_load(open('params.yaml')))"
```

---

## Getting Help

### Check Documentation

- **Main Documentation**: `docs/decontamination_optimization.md`
- **YAML Guide**: `docs/YAML_PARAMETERS_GUIDE.md`
- **Migration Guide**: `docs/MIGRATION_GUIDE.md`
- **Quick Reference**: `docs/DECONTAMINATION_QUICK_REFERENCE.md`

### Report Issues

If you encounter a bug or need help:

1. **Collect information:**
   ```bash
   # Pipeline version
   nextflow run main.nf --version
   
   # Nextflow version
   nextflow -version
   
   # Last 100 lines of log
   tail -n 100 .nextflow.log > error.log
   ```

2. **Create GitHub issue:**
   - Go to: https://github.com/gene2dis/BugBuster/issues
   - Include: error message, configuration, and logs
   - Describe: what you expected vs. what happened

3. **Include:**
   - Nextflow version
   - Pipeline version
   - Configuration (params.yaml or command)
   - Error message
   - Relevant log sections

---

## Quick Fixes Checklist

- [ ] Using latest pipeline version
- [ ] Samplesheet file exists and is formatted correctly
- [ ] FASTA files exist and are accessible
- [ ] Sufficient disk space available
- [ ] Sufficient memory allocated
- [ ] Using absolute paths for files
- [ ] YAML syntax is correct (no spaces after commas in lists)
- [ ] Work directory exists (for resume)
- [ ] Container/environment is working (Docker/Singularity)
