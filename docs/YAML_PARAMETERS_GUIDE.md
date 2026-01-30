# YAML Parameters Guide: Multiple Genomes

## Overview

This guide explains how to specify multiple contaminant genomes in YAML parameter files for BugBuster's optimized single-pass decontamination.

## Basic Syntax

### Multiple Genomes as Comma-Separated String

The recommended approach is to provide multiple FASTA files as a **comma-separated string**:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

custom_host_fasta: "/path/to/human.fasta,/path/to/mouse.fasta,/path/to/ecoli.fasta"
custom_phiX_fasta: "/path/to/phix174.fasta"
```

## Common Use Cases

### 1. Human + Mouse Samples

For samples that may contain both human and mouse contamination:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

custom_host_fasta: "/data/genomes/human_GRCh38.fasta,/data/genomes/mouse_GRCm39.fasta"
custom_phiX_fasta: "/data/genomes/phix174.fasta"
```

**Usage:**
```bash
nextflow run main.nf -params-file params.yaml
```

### 2. Multiple Host Organisms

For complex samples with multiple potential hosts:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

# Human + Mouse + Rat
custom_host_fasta: "/data/genomes/human.fasta,/data/genomes/mouse.fasta,/data/genomes/rat.fasta"
custom_phiX_fasta: "/data/genomes/phix.fasta"

# Optional: adjust resources for larger index
max_memory: "64.GB"
max_cpus: 16
```

### 3. Host + Common Lab Contaminants

Remove host plus common laboratory contaminants:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

# Human + E. coli + Yeast
custom_host_fasta: "/data/genomes/human.fasta,/data/genomes/ecoli_K12.fasta,/data/genomes/yeast_S288C.fasta"
custom_phiX_fasta: "/data/genomes/phix.fasta"
```

### 4. Host + Cloning Vectors

Remove host genome and common cloning vectors:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

# Human + pUC19 + pBR322 vectors
custom_host_fasta: "/data/genomes/human.fasta,/data/vectors/pUC19.fasta,/data/vectors/pBR322.fasta"
custom_phiX_fasta: "/data/genomes/phix.fasta"
```

### 5. Compressed FASTA Files

Mix compressed and uncompressed files (both work):

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

# Mix of .gz and uncompressed
custom_host_fasta: "/data/genomes/human.fasta.gz,/data/genomes/mouse.fasta,/data/genomes/ecoli.fasta.gz"
custom_phiX_fasta: "/data/genomes/phix174.fasta.gz"
```

## Advanced Configurations

### With Quoted Paths (Spaces in Filenames)

If your paths contain spaces, use quotes:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

# Note the single quotes around the entire string, and double quotes around each path
custom_host_fasta: '"/data/genomes/human genome.fasta","/data/genomes/mouse genome.fasta"'
```

### Pre-built Combined Index

If you've already built a combined index from multiple genomes:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

# Use pre-built index (fastest option)
custom_decontamination_index: "/data/databases/contaminants_index"
```

**How to build the index:**
```bash
# Concatenate all genomes
cat human.fasta mouse.fasta ecoli.fasta phix.fasta > contaminants.fasta

# Build index
mkdir contaminants_index
bowtie2-build --threads 16 contaminants.fasta contaminants_index/contaminants
```

### Complete Production Configuration

Full parameter set with multiple genomes:

```yaml
# Input/Output
input: "samplesheet.csv"
output: "results"

# Quality Control
quality_control: true
min_read_sample: 10000
store_clean_reads: true

# Multiple contaminant genomes
custom_host_fasta: "/data/genomes/human.fasta,/data/genomes/mouse.fasta,/data/genomes/rat.fasta"
custom_phiX_fasta: "/data/genomes/phix.fasta"

# Bowtie2 alignment parameters
bowtie_ma: 2
bowtie_mp: "6,2"
bowtie_score_min: "G,15,6"
bowtie_k: 1
bowtie_N: 1
bowtie_L: 20
bowtie_R: 2
bowtie_i: "S,1,0.75"

# Resource limits
max_memory: "128.GB"
max_cpus: 32
max_time: "240.h"

# Database storage
databases_dir: "/shared/databases"

# Assembly and profiling
assembly_mode: "assembly"
taxonomic_profiler: "sourmash"

# Feature toggles
read_arg_prediction: false
rgi_prediction: false
contig_tax_and_arg: false
include_binning: false
```

## Important Notes

### 1. Comma-Separated Format

✅ **Correct:**
```yaml
custom_host_fasta: "/path/file1.fasta,/path/file2.fasta,/path/file3.fasta"
```

❌ **Incorrect (no spaces after commas):**
```yaml
custom_host_fasta: "/path/file1.fasta, /path/file2.fasta, /path/file3.fasta"
```

### 2. YAML List Format

The pipeline currently expects comma-separated strings, **not** YAML lists:

❌ **Not supported (yet):**
```yaml
custom_host_fasta:
  - "/path/file1.fasta"
  - "/path/file2.fasta"
  - "/path/file3.fasta"
```

### 3. File Path Requirements

- Use **absolute paths** or paths relative to where you run Nextflow
- Ensure all FASTA files exist and are readable
- Both compressed (`.gz`) and uncompressed files work
- Files can be mixed (some compressed, some not)

### 4. Memory Considerations

More genomes = larger combined index = more memory needed:

```yaml
# For 2-3 genomes
max_memory: "32.GB"

# For 4-6 genomes
max_memory: "64.GB"

# For 7+ genomes
max_memory: "128.GB"
```

## Usage Examples

### Basic Usage

```bash
# Create params.yaml with your configuration
nextflow run main.nf -params-file params.yaml
```

### With Profile

```bash
# Use Docker profile
nextflow run main.nf -params-file params.yaml -profile docker

# Use Slurm profile
nextflow run main.nf -params-file params.yaml -profile slurm

# Use low disk profile
nextflow run main.nf -params-file params.yaml -profile low_disk
```

### Override Specific Parameters

```bash
# Override output directory
nextflow run main.nf -params-file params.yaml --output different_results

# Override max memory
nextflow run main.nf -params-file params.yaml --max_memory 256.GB
```

### Resume Failed Run

```bash
# Resume from last checkpoint
nextflow run main.nf -params-file params.yaml -resume
```

## Validation

### Check Your YAML Syntax

```bash
# Validate YAML syntax
python3 -c "import yaml; yaml.safe_load(open('params.yaml'))"
```

### Test File Paths

```bash
# Check if files exist
for file in $(grep "custom_.*_fasta:" params.yaml | cut -d'"' -f2 | tr ',' '\n'); do
  echo "Checking: $file"
  ls -lh "$file"
done
```

### Dry Run

```bash
# See what would be executed without running
nextflow run main.nf -params-file params.yaml -preview
```

## Troubleshooting

### Issue: "File not found" error

**Solution:** Use absolute paths or verify relative paths:
```yaml
# Absolute path (recommended)
custom_host_fasta: "/data/genomes/human.fasta,/data/genomes/mouse.fasta"

# Relative path (from where you run nextflow)
custom_host_fasta: "genomes/human.fasta,genomes/mouse.fasta"
```

### Issue: YAML parsing error

**Solution:** Check for proper quoting:
```yaml
# Correct
custom_host_fasta: "/path/file1.fasta,/path/file2.fasta"

# Also correct (with quotes)
custom_host_fasta: '/path/file1.fasta,/path/file2.fasta'
```

### Issue: High memory usage

**Solution:** Increase memory allocation:
```yaml
max_memory: "128.GB"

# Or in command line
nextflow run main.nf -params-file params.yaml --max_memory 128.GB
```

### Issue: Index building takes too long

**Solution:** Pre-build the index once:
```bash
# Build index manually
cat genome1.fasta genome2.fasta genome3.fasta > combined.fasta
bowtie2-build --threads 16 combined.fasta contaminants_index/contaminants

# Then use in YAML
custom_decontamination_index: "/path/to/contaminants_index"
```

## Complete Example Files

See the `examples/` directory for complete YAML examples:
- `params_multiple_genomes.yaml` - Multiple genome configurations
- `params_human_mouse.yaml` - Human + mouse specific example
- `params_production.yaml` - Full production configuration

## References

- **Nextflow Parameters**: https://www.nextflow.io/docs/latest/config.html#scope-params
- **YAML Syntax**: https://yaml.org/spec/1.2/spec.html
- **BugBuster Documentation**: `docs/decontamination_optimization.md`
