# YAML Parameters Guide: Multiple Genomes

## Overview

This guide explains how to specify multiple contaminant genomes (and other list-like options) in YAML parameter files for BugBuster.

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

# Larger index: the default aligner (bwa-mem3) needs ~20-25 GB RAM per Gbp
# of combined genome to build it (human + mouse + rat ~8.4 Gbp); see
# "Memory Considerations" below
max_memory: "200.GB"
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

### Paths with Spaces or Commas

Not supported: the list is split on commas, and the index build reads the staged
file names unquoted. Link such files under plain names first:

```bash
ln -s "/data/genomes/human genome.fasta" /data/genomes/human.fasta
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

# Build the bwa-mem3 index (default aligner)
mkdir contaminants_index
bwa-mem3 index -t 16 -p contaminants_index/contaminants contaminants.fasta

# Or a Bowtie2 index, used with read_aligner: "bowtie2"
mkdir bowtie2_index
bowtie2-build --threads 16 contaminants.fasta bowtie2_index/contaminants
```

### Complete Production Configuration

Full parameter set with multiple genomes:

```yaml
# Input/Output
input: "samplesheet.csv"
output: "results"

# Quality Control
quality_control: true
min_read_sample: 10000000  # optional; default 0 disables the sample filter
store_clean_reads: true

# Multiple contaminant genomes
custom_host_fasta: "/data/genomes/human.fasta,/data/genomes/mouse.fasta,/data/genomes/rat.fasta"
custom_phiX_fasta: "/data/genomes/phix.fasta"

# Aligner: bwa-mem3 (default) or bowtie2 (fallback for CPUs without AVX2
# or nodes with too little RAM for the bwa-mem3 index build)
read_aligner: "bwamem3"
decontam_min_score: 30        # bwa-mem3 only

# Bowtie2 alignment parameters (read_aligner: "bowtie2" only)
bowtie_ma: 2
bowtie_mp: "6,2"
bowtie_score_min: "G,15,6"
bowtie_k: 1
bowtie_N: 1
bowtie_L: 20
bowtie_R: 2
bowtie_i: "S,1,0.75"

# Resource limits: building the bwa-mem3 index of human + mouse + rat
# (~8.4 Gbp) requests ~200 GB (see "Memory Considerations")
max_memory: "200.GB"
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
kmarvel_prediction: false
contig_tax_and_arg: false
contig_level_functional: false
mag_level_functional: false
read_level_functional: "none"   # or "woltka", "superfocus", "humann"
superfocus_aligner: "diamond"   # SUPER-FOCUS only: or "mmseqs2"
humann_db: "v4_alpha-full"      # HUMAnN only: or "v4_alpha-ec_filtered" (smaller ChocoPhlAn)
include_binning: false
```

## Important Notes

### 1. Comma-Separated Format

Spaces around the commas are ignored, so both forms work:
```yaml
custom_host_fasta: "/path/file1.fasta,/path/file2.fasta,/path/file3.fasta"
custom_host_fasta: "/path/file1.fasta, /path/file2.fasta, /path/file3.fasta"
```

`custom_phiX_fasta` takes a single FASTA.

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

More genomes = larger combined index = more memory needed. The index is built
once and can be reused with `custom_decontamination_index`.

With the default aligner (bwa-mem3), the build needs about 20-25 GB RAM per Gbp
of combined genome (measured: 57 GB peak for GRCh38 + phiX, 3.1 Gbp). The task
requests ~24 GB per Gbp, computed from the FASTA sizes (a gzipped file counts as
4x its size), capped by `max_memory`; if that is too little, the build stops
with exit status 3. The pipeline warns at launch when `max_memory` is below
75 GB only for the built-in human host (`host_db`), not for
`custom_host_fasta`, so size `max_memory` yourself for custom genomes:

```yaml
# Human only
max_memory: "80.GB"

# Human + mouse (~5.8 Gbp)
max_memory: "150.GB"
```

With `read_aligner: "bowtie2"` the build needs far less memory (~6 GB for a
human host) but takes about an hour:

```yaml
# For 2-3 genomes
max_memory: "32.GB"

# For 4-6 genomes
max_memory: "64.GB"
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
for file in $(grep "custom_.*_fasta:" params.yaml | cut -d'"' -f2 | tr ',' '\n' | sed 's/^ *//'); do
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
# Build the bwa-mem3 index (default aligner) manually
cat genome1.fasta genome2.fasta genome3.fasta > combined.fasta
mkdir contaminants_index
bwa-mem3 index -t 16 -p contaminants_index/contaminants combined.fasta

# Then use in YAML
custom_decontamination_index: "/path/to/contaminants_index"
```

A Bowtie2 index (`bowtie2-build`) is used the same way, with
`read_aligner: "bowtie2"`; under the default aligner it is rejected at launch.
The index the pipeline builds is also published to
`<databases_dir>/bwamem3_index/contaminants_index` (bowtie2:
`<databases_dir>/bowtie_index/contaminants_index`), ready to pass here.

## Complete Example Files

See the `examples/` directory for a complete YAML example:
- `params_multiple_genomes.yaml` - Multiple genome configurations

## Binning Configuration

### Selecting Binners

Use the `binners` parameter to choose which binning tools to run. Accepts a comma-separated string with any combination of `semibin`, `metabat2`, and `comebin`. At least one is required when `include_binning: true`.

**Key rule:** Binette bin refinement is **only invoked when ≥2 binners are selected**.

#### Single Binner (Fast, Default)

```yaml
include_binning: true
binners: "semibin"
semibin_env_model: "human_gut"
```

```bash
nextflow run main.nf --input samplesheet.csv --output results \
    --include_binning true --binners semibin --semibin_env_model human_gut \
    -profile docker
```

#### Two Binners + Binette Refinement

```yaml
include_binning: true
binners: "semibin,metabat2"
semibin_env_model: "human_gut"
bin_min_completeness: 50
bin_max_contamination: 10
```

```bash
nextflow run main.nf --input samplesheet.csv --output results \
    --include_binning true --binners semibin,metabat2 \
    --semibin_env_model human_gut \
    --bin_min_completeness 50 --bin_max_contamination 10 \
    -profile docker
```

#### All Three Binners + Binette Refinement

```yaml
include_binning: true
binners: "semibin,metabat2,comebin"
semibin_env_model: "human_gut"
bin_min_completeness: 50
bin_max_contamination: 10
```

```bash
nextflow run main.nf --input samplesheet.csv --output results \
    --include_binning true --binners semibin,metabat2,comebin \
    --bin_min_completeness 50 --bin_max_contamination 10 \
    -profile docker
```

#### Complete Binning Configuration in YAML

```yaml
# Input/Output
input: "samplesheet.csv"
output: "results"

# Pipeline
assembly_mode: "assembly"
taxonomic_profiler: "sourmash"
include_binning: true

# Binner selection — ≥2 enables Binette refinement automatically
binners: "semibin,metabat2"

# SemiBin
semibin_env_model: "human_gut"

# Refinement thresholds (only used when ≥2 binners selected)
bin_min_completeness: 50
bin_max_contamination: 10

# MetaBAT2 contig filter (only used when metabat2 selected)
metabat_minContig: 2500
```

### Valid Binner Values

| Value | Tool | Notes |
|-------|------|-------|
| `semibin` | SemiBin2 | Fast; requires `--semibin_env_model` |
| `metabat2` | MetaBAT2 | Requires contig depth calculation |
| `comebin` | COMEBin | Slow on CPU; best with GPU |

### Important Notes

- `binners` accepts a **comma-separated string**: `"semibin,metabat2"`; spaces around the commas are ignored in YAML (on the command line, write it without spaces)
- `bin_min_completeness` and `bin_max_contamination` are **ignored** when only one binner is selected (`metawrap_completeness` / `metawrap_contamination` are deprecated aliases)
- CheckM2 and GTDB-Tk **always run** regardless of binner selection (with ≥2 binners CheckM2 runs inside Binette)

---

## References

- **Nextflow Parameters**: https://www.nextflow.io/docs/latest/config.html#scope-params
- **YAML Syntax**: https://yaml.org/spec/1.2/spec.html
- **BugBuster Documentation**: [`manual.md`](manual.md), [`parameters.md`](parameters.md)
