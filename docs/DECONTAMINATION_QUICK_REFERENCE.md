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
The default, bwa-mem3, builds its index for a human host in a few minutes with ~75 GB RAM; mapping against it needs ~14 GB per task. A pair counts as clean only when **both** mates are unmapped; each mate must reach `--decontam_min_score` on its own (bwa-mem3 runs with pairing off, `-P`). `--custom_decontamination_index` must then point to a bwa-mem3 or bwa-mem2 index (`*.bwt.2bit.64`); a Bowtie2 index is rejected at launch.

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

# Build the bwa-mem3 index (default aligner; ~75 GB RAM for a human host)
mkdir contaminants_index
bwa-mem3 index -t 16 -p contaminants_index/contaminants contaminants.fasta
nextflow run main.nf --custom_decontamination_index contaminants_index

# Or a Bowtie2 index for the fallback
bowtie2-build contaminants.fasta bowtie2_index/contaminants
nextflow run main.nf --read_aligner bowtie2 --custom_decontamination_index bowtie2_index
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
<databases_dir>/bwamem3_index/contaminants_index/  # default: <output>/../databases
<databases_dir>/bowtie_index/contaminants_index/   # under --read_aligner bowtie2
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
The bwa-mem3 index build for a human host needs ~75 GB RAM. With a lower
`--max_memory` the pipeline warns at launch and the build stops with exit status 3.
Pass a pre-built index, or use the bowtie2 fallback:
```bash
nextflow run main.nf --max_memory 64.GB --read_aligner bowtie2
```

### "Illegal instruction" in a BWAMEM3_* task
The bioconda bwa-mem3 build needs a CPU with AVX2 (x86-64-v3). Use `--read_aligner bowtie2`.

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
