# Decontamination Quick Reference

## TL;DR

BugBuster removes phiX and host contamination in a **single pass** against one combined bwa-mem3 index, built automatically from the phiX and host genome FASTAs. Building it for a human host needs ~75 GB RAM; `--read_aligner bowtie2` is the fallback for smaller nodes and for CPUs without AVX2.

---

## Quick Start

### Basic Usage
```bash
nextflow run main.nf --input samples.csv --output results -profile docker
```

### With Pre-built Index (Fastest)
```bash
nextflow run main.nf \
  --input samples.csv \
  --output results \
  --custom_decontamination_index /path/to/contaminants_index \
  -profile docker
```

### Custom Host Genome
```bash
nextflow run main.nf \
  --input samples.csv \
  --output results \
  --custom_host_fasta /path/to/host.fasta \
  -profile docker
```

### Aligners
The default, bwa-mem3, builds its index for a human host in a few minutes with ~75 GB RAM; mapping against it peaked at ~16 GB per task in the [benchmark](benchmarks/bwamem3_vs_bowtie2.md) (the task requests 72 GB, capped by `--max_memory`). A pair counts as clean only when **both** mates are unmapped; each mate must reach `--decontam_min_score` on its own (bwa-mem3 runs with pairing off, `-P`). `--custom_decontamination_index` must then point to a bwa-mem3 or bwa-mem2 index (`*.bwt.2bit.64`); a Bowtie2 index is rejected at launch.

### With bowtie2 (fallback: `--read_aligner bowtie2`)
```bash
nextflow run main.nf \
  --input samples.csv \
  --output results \
  --read_aligner bowtie2 \
  -profile docker
```
Use it on CPUs without AVX2, on nodes with less than ~75 GB RAM, or to keep using a Bowtie2 index from an earlier release (`--custom_decontamination_index`). Building the human Bowtie2 index takes about an hour with ~6 GB RAM. bowtie2 keeps a pair when it has no concordant alignment (`--un-conc`), so a pair with only one host mate survives, as do host pairs aligned with an insert above 500 bp. See [`--read_aligner`](parameters.md#--read_aligner).

---

## Key Parameters

| Parameter | Purpose | Example |
|-----------|---------|---------|
| `--custom_decontamination_index` | Pre-built combined index | `/data/contaminants_index` |
| `--custom_phiX_fasta` | Custom phiX FASTA | `/data/phix.fasta` |
| `--custom_host_fasta` | Custom host FASTA, or a comma-separated list | `/data/host.fasta` |
| `--quality_control` | Enable/disable QC | `true` (default) |
| `--store_clean_reads` | Store clean reads permanently | `false` (default) |

---

## Common Workflows

### 1. Default (Human + PhiX)
```bash
nextflow run main.nf --input samples.csv --output results -profile docker
```
Automatically downloads and uses human CHM13 + phiX174.

### 2. Non-Human Host (e.g. Mouse)
```bash
nextflow run main.nf \
  --input samples.csv \
  --output results \
  --custom_host_fasta /path/to/mouse_genome.fasta.gz \
  -profile docker
```
(`--host_db` has a single built-in entry, `human`; other hosts are supplied as FASTA.)

### 3. Multiple Contaminants
```bash
nextflow run main.nf \
  --input samples.csv \
  --output results \
  --custom_host_fasta "human.fasta,mouse.fasta,ecoli.fasta" \
  -profile docker
```
All the FASTAs go into the one combined index. The bwa-mem3 build requests ~24 GB RAM
per Gbp of combined genome, so raise `--max_memory` for several large genomes (see
[Memory Considerations](#4-memory-considerations) below).

### 4. Low Disk Space
Deletes the work directory after a successful run and keeps the clean reads in the output
dir. It does not lower peak usage during the run, and a completed run cannot be resumed
(see [`DISK_OPTIMIZATION.md`](DISK_OPTIMIZATION.md)).
```bash
nextflow run main.nf \
  --input samples.csv \
  --output results \
  -profile docker,low_disk
```

---

## Building Combined Index

The one recipe for building a combined (host + phiX) index by hand; the pipeline builds
and stores it itself when given the FASTA files, so this is for an index shared between
`--databases_dir`s or built on another machine.

```bash
# Concatenate FASTA files
cat phix.fasta host.fasta > contaminants.fasta

# Build the bwa-mem3 index (default aligner; ~75 GB RAM for a human host)
mkdir contaminants_index
bwa-mem3 index -t 16 -p contaminants_index/contaminants contaminants.fasta
nextflow run main.nf --input samples.csv --output results -profile docker \
    --custom_decontamination_index contaminants_index

# Or a Bowtie2 index for the fallback
mkdir bowtie2_index
bowtie2-build contaminants.fasta bowtie2_index/contaminants
nextflow run main.nf --input samples.csv --output results -profile docker \
    --read_aligner bowtie2 --custom_decontamination_index bowtie2_index
```

---

## Multiple Host Genomes in a Params File

Several contaminant genomes (and other list-like options) in a `-params-file` YAML
(general rules for params files: [`YAML_PARAMETERS_GUIDE.md`](YAML_PARAMETERS_GUIDE.md)).
Run each example with `nextflow run main.nf -params-file params.yaml -profile docker`.

### Basic Syntax

#### Multiple Genomes as Comma-Separated String

The recommended approach is to provide multiple FASTA files as a **comma-separated string**:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

custom_host_fasta: "/path/to/human.fasta,/path/to/mouse.fasta,/path/to/ecoli.fasta"
custom_phiX_fasta: "/path/to/phix174.fasta"
```

### Common Use Cases

#### 1. Human + Mouse Samples

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
nextflow run main.nf -params-file params.yaml -profile docker
```

#### 2. Multiple Host Organisms

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

#### 3. Host + Common Lab Contaminants

Remove host plus common laboratory contaminants:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

# Human + E. coli + Yeast
custom_host_fasta: "/data/genomes/human.fasta,/data/genomes/ecoli_K12.fasta,/data/genomes/yeast_S288C.fasta"
custom_phiX_fasta: "/data/genomes/phix.fasta"
```

#### 4. Host + Cloning Vectors

Remove host genome and common cloning vectors:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

# Human + pUC19 + pBR322 vectors
custom_host_fasta: "/data/genomes/human.fasta,/data/vectors/pUC19.fasta,/data/vectors/pBR322.fasta"
custom_phiX_fasta: "/data/genomes/phix.fasta"
```

#### 5. Compressed FASTA Files

Mix compressed and uncompressed files (both work):

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

# Mix of .gz and uncompressed
custom_host_fasta: "/data/genomes/human.fasta.gz,/data/genomes/mouse.fasta,/data/genomes/ecoli.fasta.gz"
custom_phiX_fasta: "/data/genomes/phix174.fasta.gz"
```

### Advanced Configurations

#### Paths with Spaces or Commas

Not supported: the list is split on commas, and the index build reads the staged
file names unquoted. Link such files under plain names first:

```bash
ln -s "/data/genomes/human genome.fasta" /data/genomes/human.fasta
```

#### Pre-built Combined Index

If you've already built a combined index from multiple genomes:

```yaml
input: "samplesheet.csv"
output: "results"
quality_control: true

# Use pre-built index (fastest option)
custom_decontamination_index: "/data/databases/contaminants_index"
```

How to build it: [Building Combined Index](#building-combined-index) above (concatenate
all the genomes and phiX into one FASTA first).

#### Complete Production Configuration

Full parameter set with multiple genomes:

```yaml
# Input/Output
input: "samplesheet.csv"
output: "results"

# Quality Control
quality_control: true
min_read_sample: 20000000  # optional; single reads (R1 + R2) after fastp, so 10 M pairs; default 0 disables the filter
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
taxonomic_profiler: "sylph"     # default; or "kraken2", "sourmash", "none"
nonpareil: true                 # default: sequencing coverage per sample

# Feature toggles
args_oap_prediction: false
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

### Important Notes

#### 1. Comma-Separated Format

Spaces around the commas are ignored, so both forms work:
```yaml
custom_host_fasta: "/path/file1.fasta,/path/file2.fasta,/path/file3.fasta"
custom_host_fasta: "/path/file1.fasta, /path/file2.fasta, /path/file3.fasta"
```

`custom_phiX_fasta` takes a single FASTA.

#### 2. YAML List Format

The pipeline currently expects comma-separated strings, **not** YAML lists:

❌ **Not supported (yet):**
```yaml
custom_host_fasta:
  - "/path/file1.fasta"
  - "/path/file2.fasta"
  - "/path/file3.fasta"
```

#### 3. File Path Requirements

- Use **absolute paths** or paths relative to where you run Nextflow
- Ensure all FASTA files exist and are readable
- Both compressed (`.gz`) and uncompressed files work
- Files can be mixed (some compressed, some not)

#### 4. Memory Considerations

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

## Output Files

### Clean Reads (with `--store_clean_reads`)
```
results/clean_reads/sample1/
├── sample1_R1_clean.fastq.gz
├── sample1_R2_clean.fastq.gz
└── sample1_Singleton_clean.fastq.gz (if present)
```

### Reports
```
results/01_quality_control/summary/Reads_report.csv
```

### Combined Index (Stored)
```
<databases_dir>/bwamem3_index/phiX174_human/contaminants_index/  # default: <output>/../databases; reused by later runs
<databases_dir>/bowtie_index/phiX174_human/contaminants_index/   # under --read_aligner bowtie2
<databases_dir>/bwamem3_index/contaminants_index/                # copy of an index built from --custom_host_fasta / --custom_phiX_fasta
```

---

## Performance Tips

1. **Reuse the index**: keep the same `--databases_dir` and later runs skip the build; an index built from `--custom_host_fasta` is rebuilt every run unless you pass its copy (`<databases_dir>/bwamem3_index/contaminants_index`) with `--custom_decontamination_index`
2. **Use `--store_clean_reads true`** to keep clean reads in the output dir (`-resume` remains the task cache)
3. **Use `-profile low_disk`** to free the work directory after a successful run (it does not lower peak usage; put `-work-dir` on a large filesystem)

---

## Troubleshooting

Decontamination problems are in [`troubleshooting.md`](troubleshooting.md#decontamination-issues):
the index build stopping for memory (exit status 3), "Illegal instruction" (exit status
132, a CPU without AVX2), rejected indexes, multiple host genomes and slow index builds.
To check how many reads each step removed, look at
`results/01_quality_control/summary/Reads_report.csv`.

---

## How It Works

```
Reads → Combined phiX + Host Decontamination (one bwa-mem3 or Bowtie2 pass) → Clean Reads
```

A single alignment against the combined index replaces separate host and phiX
passes: fewer temp files, less disk, one decontamination report per sample.

---

## Documentation

- **Manual**: [`manual.md`](manual.md)
- **Parameter reference**: [`parameters.md`](parameters.md)
- **Troubleshooting**: [`troubleshooting.md`](troubleshooting.md)
- **Examples**: `examples/decontamination_examples.sh`

---

## Support

- GitHub: [https://github.com/gene2dis/BugBuster](https://github.com/gene2dis/BugBuster)
- Issues: [https://github.com/gene2dis/BugBuster/issues](https://github.com/gene2dis/BugBuster/issues)
