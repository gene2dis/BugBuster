# Decontamination Quick Reference

## TL;DR

BugBuster removes phiX and host contamination in a **single pass** against one combined Bowtie2 index, built automatically from the phiX and host genome FASTAs.

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

---

## Key Parameters

| Parameter | Purpose | Example |
|-----------|---------|---------|
| `--custom_decontamination_index` | Pre-built combined index | `/data/contaminants_index` |
| `--custom_phiX_fasta` | Custom phiX FASTA | `/data/phix.fasta` |
| `--custom_host_fasta` | Custom host FASTA(s) | `/data/host.fasta` |
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

```bash
# Concatenate FASTA files
cat phix.fasta host.fasta > contaminants.fasta

# Build index
mkdir contaminants_index
bowtie2-build contaminants.fasta contaminants_index/contaminants

# Use in pipeline
nextflow run main.nf --custom_decontamination_index contaminants_index
```

---

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

### Combined Index (Cached)
```
<databases_dir>/bowtie_index/contaminants_index/   # default: <output>/../databases
```

---

## Performance Tips

1. **Pre-build index** for multiple runs
2. **Use `--store_clean_reads true`** to keep clean reads in the output dir (`-resume` remains the task cache)
3. **Use `-profile low_disk`** to free the work directory after a successful run (it does not lower peak usage; put `-work-dir` on a large filesystem)
4. **Share `--databases_dir`** across projects

---

## Troubleshooting

### High Memory Usage
```bash
nextflow run main.nf --max_memory 64.GB
```

### Resume Failed Run
```bash
nextflow run main.nf -resume
```

### Check Statistics
```bash
cat results/01_quality_control/summary/Reads_report.csv
```

---

## How It Works

```
Reads → Combined phiX + Host Decontamination (one Bowtie2 pass) → Clean Reads
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
