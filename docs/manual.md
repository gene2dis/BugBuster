# BugBuster User Manual

**Bacterial Unraveling and metaGenomic Binning with Up-Scale Throughput, Efficient and Reproducible**

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Requirements](#2-requirements)
3. [Installation](#3-installation)
4. [Database Management](#4-database-management)
5. [Input Preparation](#5-input-preparation)
6. [Pipeline Parameters](#6-pipeline-parameters)
7. [Running the Pipeline](#7-running-the-pipeline)
8. [Output Structure](#8-output-structure)
9. [Usage Examples](#9-usage-examples)
10. [Advanced Configuration](#10-advanced-configuration)
11. [Troubleshooting](#11-troubleshooting)

---

## 1. Introduction

BugBuster is a comprehensive Nextflow pipeline for microbial metagenomic analysis and antimicrobial resistance gene (ARG) prediction. Built following bioinformatics best practices, it provides:

- **Read Quality Control**: Filtering, trimming, and host decontamination
- **Taxonomic Profiling**: Species identification using Kraken2 or Sourmash
- **Genome Assembly**: Per-sample or co-assembly using MEGAHIT
- **Metagenomic Binning**: MAG recovery with MetaBAT2, SemiBin, and COMEBin
- **Bin Refinement**: Quality improvement with MetaWRAP
- **Bin Quality Assessment**: Completeness and contamination with CheckM2
- **Taxonomic Classification**: Bin taxonomy with GTDB-TK
- **ARG Prediction**: At read, contig, and bin levels
- **Functional Annotation**: With MetaCerberus

### Pipeline Workflow

```
Reads → QC (FastP) → Host Removal (Bowtie2) → Taxonomy (Kraken2/Sourmash)
                                            ↓
                                      Assembly (MEGAHIT)
                                            ↓
                                      Contig Filtering (BBMap)
                                            ↓
                              ┌─────────────┴─────────────┐
                              ↓                           ↓
                      Binning (3 tools)           Contig Analysis
                              ↓                    (BlobTools, DeepARG)
                      Refinement (MetaWRAP)
                              ↓
                      Quality (CheckM2)
                              ↓
                      Taxonomy (GTDB-TK)
```

---

## 2. Requirements

### Software Requirements

| Software | Minimum Version | Purpose |
|----------|-----------------|---------|
| Nextflow | ≥23.04.0 | Workflow engine |
| Java | 11-17 | Nextflow runtime |
| Container runtime | - | Docker, Singularity, or Podman |

### Hardware Requirements

| Analysis Type | Memory | CPUs | Storage |
|---------------|--------|------|---------|
| QC + Taxonomy only | 16 GB | 4 | 50 GB |
| QC + Taxonomy + Assembly | 64 GB | 16 | 100 GB |
| Full pipeline with binning | 128 GB | 16 | 200 GB |
| Large datasets / co-assembly | 256 GB | 32 | 500 GB |

---

## 3. Installation

### Install Nextflow

```bash
# Using curl
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

### Install Container Runtime

Choose one of the following:

**Docker:**
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install docker.io
sudo usermod -aG docker $USER

# Verify
docker run hello-world
```

**Singularity:**
```bash
# Ubuntu/Debian
sudo apt-get install singularity-container

# Verify
singularity --version
```

### Clone the Pipeline

```bash
git clone https://github.com/gene2dis/BugBuster.git
cd BugBuster
```

---

## 4. Database Management

BugBuster automatically downloads required databases on first use. Databases are stored in `<output_dir>/downloaded_db/` with symbolic links for reuse.

### Automatic Download Databases

| Database | Size | Used For | Download Trigger |
|----------|------|----------|------------------|
| **phiX174 Index** | 8.1 MB | PhiX contamination removal | `quality_control=true` |
| **Human Host Index** | 4.1 GB | Host read removal | `quality_control=true` |
| **Kraken2 Standard-8** | 7.5 GB | Taxonomic profiling | `taxonomic_profiler='kraken2'` |
| **Kraken2 GTDB r220** | 497 GB | Taxonomic profiling | `kraken2_db='gtdb_220'` |
| **Sourmash GTDB r220** | 17 GB | Taxonomic profiling | `taxonomic_profiler='sourmash'` |
| **NCBI Taxdump** | 448 MB | BlobTools taxonomy | `contig_tax_and_arg=true` |
| **NCBI NT** | 434 GB | Contig BLAST | `contig_tax_and_arg=true` |
| **DeepARG** | 4.8 GB | Contig ARG prediction | `contig_tax_and_arg=true` |
| **KARGA (MEGARes)** | 9.2 MB | Read ARG prediction | `read_arg_prediction=true` |
| **KARGVA** | 1.5 MB | Read ARG variant prediction | `read_arg_prediction=true` |
| **CheckM2** | 2.9 GB | Bin quality assessment | `include_binning=true` |
| **GTDB-TK r220** | 109 GB | Bin taxonomic classification | `include_binning=true` |

### Manual Database Download

For shared HPC systems or repeated runs, pre-download databases to shared storage:

```bash
# Create database directory
mkdir -p /shared/databases/bugbuster

# Download Kraken2 Standard-8
wget -O - https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz | \
    tar -xzf - -C /shared/databases/bugbuster/kraken2_standard8

# Download Sourmash GTDB
wget -P /shared/databases/bugbuster/sourmash/ \
    https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs220/gtdb-reps-rs220-k31.zip
wget -P /shared/databases/bugbuster/sourmash/ \
    https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs220/gtdb-rs220.lineages.reps.csv

# Download CheckM2
wget -O /shared/databases/bugbuster/checkm2_db.tar.gz \
    "https://zenodo.org/records/14897628/files/checkm2_database.tar.gz?download=1"
tar -xzf /shared/databases/bugbuster/checkm2_db.tar.gz -C /shared/databases/bugbuster/

# Download GTDB-TK (large download)
wget -O /shared/databases/bugbuster/gtdbtk_r220.tar.gz \
    https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz
tar -xzf /shared/databases/bugbuster/gtdbtk_r220.tar.gz -C /shared/databases/bugbuster/

# Download human host Bowtie2 index
wget -O /shared/databases/bugbuster/chm13_plusY.zip \
    https://genome-idx.s3.amazonaws.com/bt/chm13.draft_v1.0_plusY.zip
unzip /shared/databases/bugbuster/chm13_plusY.zip -d /shared/databases/bugbuster/host_index/
```

### Using Custom Databases

Specify custom database paths to skip automatic downloads:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --custom_kraken_db /shared/databases/bugbuster/kraken2_standard8 \
    --custom_bowtie_host_index /shared/databases/bugbuster/host_index \
    --custom_checkm2_db /shared/databases/bugbuster/checkm2/uniref100.KO.1.dmnd \
    --custom_gtdbtk_db /shared/databases/bugbuster/gtdbtk_r220 \
    -profile docker
```

### Using Custom Host Genomes

To filter reads from non-human hosts, build a custom Bowtie2 index:

```bash
# Download your host genome
wget -O host_genome.fasta.gz <URL_TO_HOST_GENOME>
gunzip host_genome.fasta.gz

# Build Bowtie2 index
bowtie2-build host_genome.fasta host_index/host

# Use in pipeline
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --custom_bowtie_host_index ./host_index \
    -profile docker
```

---

## 5. Input Preparation

### Samplesheet Format

Create a CSV file with your sample information:

```csv
sample,r1,r2,s
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,
sample3,/path/to/sample3_R1.fastq.gz,/path/to/sample3_R2.fastq.gz,/path/to/sample3_singletons.fastq.gz
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier (alphanumeric, underscores allowed) |
| `r1` | Yes | Absolute path to forward reads (R1) in FASTQ/FASTQ.GZ format |
| `r2` | Yes | Absolute path to reverse reads (R2) in FASTQ/FASTQ.GZ format |
| `s` | No | Absolute path to singleton reads (optional, leave empty if none) |

### Input Requirements

- **Format**: FASTQ or compressed FASTQ (.gz)
- **Paired-end**: Required (R1 and R2)
- **Naming**: Sample names must be unique
- **Paths**: Use absolute paths for cloud/HPC execution

### Cloud Storage Paths

For cloud execution, use appropriate URI schemes:

```csv
sample,r1,r2,s
sample1,s3://bucket/data/sample1_R1.fastq.gz,s3://bucket/data/sample1_R2.fastq.gz,
sample2,gs://bucket/data/sample2_R1.fastq.gz,gs://bucket/data/sample2_R2.fastq.gz,
sample3,az://container/data/sample3_R1.fastq.gz,az://container/data/sample3_R2.fastq.gz,
```

---

## 6. Pipeline Parameters

### 6.1 Input/Output Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | *required* | Path to samplesheet CSV file |
| `--output` | *required* | Output directory for results |
| `--publish_dir_mode` | `copy` | How to save results: `copy`, `symlink`, `link`, `move` |

### 6.2 Pipeline Execution Options

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `--quality_control` | `true` | `true`, `false` | Enable QC and host filtering |
| `--assembly_mode` | `assembly` | `assembly`, `coassembly`, `none` | Genome assembly strategy |
| `--taxonomic_profiler` | `sourmash` | `kraken2`, `sourmash`, `none` | Taxonomic profiling tool |
| `--include_binning` | `false` | `true`, `false` | Enable binning and refinement |
| `--read_arg_prediction` | `false` | `true`, `false` | Read-level ARG prediction |
| `--contig_tax_and_arg` | `false` | `true`, `false` | Contig taxonomy and ARG prediction |
| `--contig_level_metacerberus` | `false` | `true`, `false` | Functional annotation with MetaCerberus |
| `--arg_bin_clustering` | `false` | `true`, `false` | ARG clustering for HGT inference |
| `--min_read_sample` | `0` | Integer ≥ 0 | Minimum reads required after QC |

### 6.3 Database Selection Options

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `--phiX_index` | `phiX174` | `phiX174` | PhiX genome for decontamination |
| `--host_db` | `human` | `human` | Host genome for read removal |
| `--kraken2_db` | `standard-8` | `standard-8`, `gtdb_220` | Kraken2 database version |
| `--sourmash_db` | `gtdb_220_k31` | `gtdb_220_k31` | Sourmash database version |
| `--checkm2_db` | `v3` | `v3` | CheckM2 database version |
| `--gtdbtk_db` | `release_220` | `release_220` | GTDB-TK database release |

### 6.4 Custom Database Paths

| Parameter | Description |
|-----------|-------------|
| `--custom_phiX_index` | Path to custom PhiX Bowtie2 index directory |
| `--custom_bowtie_host_index` | Path to custom host Bowtie2 index directory |
| `--custom_kraken_db` | Path to custom Kraken2 database directory |
| `--custom_sourmash_db` | List of paths: `["kmer.zip", "lineages.csv"]` |
| `--custom_checkm2_db` | Path to CheckM2 database file |
| `--custom_gtdbtk_db` | Path to GTDB-TK database directory |
| `--custom_deeparg_db` | Path to DeepARG database directory |
| `--custom_blast_db` | Path to BLAST NT database directory |
| `--custom_taxdump_files` | Path to NCBI taxdump directory |
| `--custom_karga_db` | Path to KARGA database FASTA file |
| `--custom_kargva_db` | Path to KARGVA database FASTA file |

### 6.5 FastP Quality Filtering Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastp_n_base_limit` | `5` | Maximum N bases allowed per read |
| `--fastp_unqualified_percent_limit` | `10` | Max percentage of unqualified bases |
| `--fastp_qualified_quality_phred` | `20` | Phred quality threshold for qualified bases |
| `--fastp_cut_front_window_size` | `4` | Sliding window size for front trimming |
| `--fastp_cut_front_mean_quality` | `20` | Mean quality threshold for front trimming |
| `--fastp_cut_right_window_size` | `4` | Sliding window size for right trimming |
| `--fastp_cut_right_mean_quality` | `20` | Mean quality threshold for right trimming |

### 6.6 Taxonomy Profiling Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--kraken_confidence` | `0.1` | Kraken2 confidence threshold (0-1) |
| `--kraken_db_used` | `gtdb_release207` | Database name for reports |
| `--bracken_read_len` | `150` | Read length for Bracken estimation |
| `--bracken_tax_level` | `S` | Taxonomic level: D, P, C, O, F, G, S |
| `--sourmash_db_name` | `gtdb_release_220` | Database name for Sourmash reports |
| `--sourmash_tax_rank` | `species` | Rank: `genus`, `species`, `strain` |

### 6.7 Assembly Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bbmap_lenght` | `1000` | Minimum contig length after filtering |

### 6.8 Binning Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metabat_minContig` | `2500` | Minimum contig length for MetaBAT2 |
| `--metawrap_completeness` | `50` | Minimum bin completeness (%) |
| `--metawrap_contamination` | `10` | Maximum bin contamination (%) |
| `--semibin_env_model` | `human_gut` | SemiBin environment model |

**SemiBin Environment Models:**
- `human_gut`, `dog_gut`, `cat_gut`, `mouse_gut`, `pig_gut`
- `human_oral`, `chicken_caecum`
- `ocean`, `soil`, `wastewater`, `built_environment`
- `global` (for mixed/unknown environments)

### 6.9 DeepARG Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--deeparg_min_prob` | `0.8` | Minimum probability threshold |
| `--deeparg_arg_alignment_identity` | `50` | Minimum alignment identity (%) |
| `--deeparg_arg_alignment_evalue` | `1e-10` | Maximum E-value |
| `--deeparg_arg_alignment_overlap` | `0.8` | Minimum alignment overlap |
| `--deeparg_model_version` | `v2` | DeepARG model version |

### 6.10 MetaCerberus Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metacerberus_hmm` | `"KOFam_all, COG, VOG, PHROG, CAZy"` | HMM databases to use |
| `--metacerberus_minscore` | `25` | Minimum HMM score |
| `--metacerberus_evalue` | `1e-09` | Maximum E-value |

**Available HMM Databases:** `KOFam_all`, `KOFam_eukaryote`, `KOFam_prokaryote`, `COG`, `VOG`, `PHROG`, `CAZy`

### 6.11 Resource Limit Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `240.h` | Maximum time per process |

---

## 7. Running the Pipeline

### 7.1 Basic Execution

```bash
# With Docker
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker

# With Singularity
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile singularity
```

### 7.2 Available Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers |
| `podman` | Run with Podman containers |
| `apptainer` | Run with Apptainer containers |
| `conda` | Run with Conda environments |
| `slurm` | Submit jobs to SLURM scheduler |
| `slurm_singularity` | SLURM with Singularity |
| `aws` | Run on AWS Batch |
| `gcp` | Run on Google Cloud |
| `azure` | Run on Azure Batch |
| `test` | Run with minimal test data |

**Combine profiles as needed:**
```bash
-profile slurm,singularity
-profile aws,docker
```

### 7.3 Resume Failed Runs

Nextflow caches completed tasks. Resume from the last successful step:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker \
    -resume
```

### 7.4 Background Execution

For long-running analyses:

```bash
# Run in background
nohup nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker \
    > pipeline.log 2>&1 &

# Or use screen/tmux
screen -S bugbuster
nextflow run main.nf --input samplesheet.csv --output ./results -profile docker
# Ctrl+A, D to detach
# screen -r bugbuster to reattach
```

### 7.5 Display Help

```bash
nextflow run main.nf --help
```

---

## 8. Output Structure

```
results/
├── pipeline_info/                    # Execution reports
│   ├── execution_report_*.html       # Resource usage report
│   ├── execution_timeline_*.html     # Timeline visualization
│   ├── execution_trace_*.txt         # Task trace log
│   └── pipeline_dag_*.svg            # Pipeline DAG
├── downloaded_db/                    # Downloaded databases (symlinks)
├── qc/                               # Quality control results
│   ├── fastp/                        # FastP reports per sample
│   └── reads_report.html             # Aggregated read statistics
├── taxonomy/                         # Taxonomic profiling
│   ├── kraken2/ or sourmash/         # Per-sample results
│   ├── phyloseq_object.rds           # R Phyloseq object
│   └── tax_report.html               # Taxonomy summary report
├── assembly/                         # Genome assembly
│   ├── {sample}/                     # Per-sample assemblies
│   │   ├── contigs.fa                # Filtered contigs
│   │   └── megahit/                  # MEGAHIT output
│   └── coassembly/                   # Co-assembly (if enabled)
├── binning/                          # Metagenomic bins
│   ├── {sample}/
│   │   ├── metabat2/                 # MetaBAT2 bins
│   │   ├── semibin/                  # SemiBin bins
│   │   ├── comebin/                  # COMEBin bins
│   │   └── metawrap/                 # Refined bins
│   └── reports/
│       ├── bin_quality_report.html   # CheckM2 quality summary
│       └── bin_taxonomy_report.html  # GTDB-TK taxonomy summary
├── args/                             # ARG predictions
│   ├── read_level/                   # KARGA/KARGVA results
│   │   ├── arg_norm_report.html      # Normalized ARG counts
│   │   └── args_oap/                 # ARGs-OAP results
│   └── contig_level/                 # DeepARG results
│       ├── deeparg/                  # Per-sample predictions
│       └── arg_contig_report.html    # Summary report
├── blobtools/                        # Contig taxonomy (BlobTools)
│   ├── {sample}/                     # Per-sample blob tables
│   └── blobplot.html                 # Aggregated visualization
└── functional/                       # Functional annotation
    └── metacerberus/                 # MetaCerberus results
```

---

## 9. Usage Examples

### Example 1: Quick Taxonomy Profiling

Minimal analysis with QC and taxonomy only:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --assembly_mode none \
    --taxonomic_profiler sourmash \
    -profile docker
```

### Example 2: Standard Metagenomic Analysis

QC, taxonomy, and per-sample assembly:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --quality_control true \
    --assembly_mode assembly \
    --taxonomic_profiler kraken2 \
    --kraken2_db standard-8 \
    -profile docker
```

### Example 3: Full Pipeline with Binning

Complete analysis including MAG recovery:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --quality_control true \
    --assembly_mode assembly \
    --taxonomic_profiler sourmash \
    --include_binning true \
    --semibin_env_model human_gut \
    --metawrap_completeness 50 \
    --metawrap_contamination 10 \
    -profile singularity
```

### Example 4: ARG-Focused Analysis

Comprehensive ARG prediction at all levels:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --quality_control true \
    --assembly_mode assembly \
    --taxonomic_profiler sourmash \
    --read_arg_prediction true \
    --contig_tax_and_arg true \
    --include_binning true \
    --arg_bin_clustering true \
    -profile docker
```

### Example 5: Co-Assembly for Related Samples

Pool samples for improved assembly:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --assembly_mode coassembly \
    --include_binning true \
    --taxonomic_profiler kraken2 \
    --kraken2_db gtdb_220 \
    -profile singularity
```

### Example 6: HPC Cluster Execution (SLURM)

```bash
nextflow run main.nf \
    --input /scratch/user/samplesheet.csv \
    --output /scratch/user/results \
    --quality_control true \
    --assembly_mode assembly \
    --include_binning true \
    --max_cpus 32 \
    --max_memory 256.GB \
    -profile slurm,singularity
```

### Example 7: AWS Batch Execution

```bash
export AWS_ACCESS_KEY_ID="your_key"
export AWS_SECRET_ACCESS_KEY="your_secret"

nextflow run main.nf \
    --input s3://bucket/samplesheet.csv \
    --output s3://bucket/results \
    --quality_control true \
    --assembly_mode assembly \
    --taxonomic_profiler sourmash \
    -profile aws,docker \
    -work-dir s3://bucket/work
```

### Example 8: Using Pre-downloaded Databases

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --quality_control true \
    --assembly_mode assembly \
    --taxonomic_profiler kraken2 \
    --include_binning true \
    --custom_bowtie_host_index /shared/db/host_index \
    --custom_kraken_db /shared/db/kraken2_standard8 \
    --custom_checkm2_db /shared/db/checkm2/uniref100.KO.1.dmnd \
    --custom_gtdbtk_db /shared/db/gtdbtk_r220 \
    -profile singularity
```

### Example 9: Non-Human Host Analysis

For mouse gut microbiome:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --custom_bowtie_host_index /path/to/mouse_bowtie2_index \
    --semibin_env_model mouse_gut \
    --include_binning true \
    -profile docker
```

### Example 10: Environmental Samples (Ocean)

```bash
nextflow run main.nf \
    --input ocean_samples.csv \
    --output ./results \
    --quality_control true \
    --assembly_mode assembly \
    --taxonomic_profiler sourmash \
    --include_binning true \
    --semibin_env_model ocean \
    --custom_bowtie_host_index null \
    -profile singularity
```

---

## 10. Advanced Configuration

### 10.1 Custom Configuration File

Create a custom config file for your environment:

```groovy
// my_config.config
params {
    max_cpus   = 32
    max_memory = '256.GB'
    max_time   = '120.h'
    
    // Pre-downloaded databases
    custom_kraken_db         = '/shared/db/kraken2_standard8'
    custom_checkm2_db        = '/shared/db/checkm2/uniref100.KO.1.dmnd'
    custom_gtdbtk_db         = '/shared/db/gtdbtk_r220'
    custom_bowtie_host_index = '/shared/db/host_index'
}

singularity {
    cacheDir = '/shared/singularity_cache'
}
```

Run with custom config:
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -c my_config.config \
    -profile singularity
```

### 10.2 Institutional Profile

For shared HPC systems, create a reusable institutional profile:

```groovy
// conf/my_institution.config
params {
    max_cpus   = 48
    max_memory = '512.GB'
}

process {
    executor = 'slurm'
    queue    = 'general'
    
    withLabel:process_high {
        queue = 'highmem'
    }
}

singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = '/shared/containers'
}
```

Add to `nextflow.config`:
```groovy
profiles {
    my_institution {
        includeConfig 'conf/my_institution.config'
    }
}
```

### 10.3 Seqera Platform Integration

Monitor runs with Seqera Platform:

```bash
export TOWER_ACCESS_TOKEN="your_token"

nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker \
    -with-tower
```

---

## 11. Troubleshooting

### Common Issues

#### Out of Memory Errors

```bash
# Increase memory limits
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --max_memory 256.GB \
    -profile docker
```

#### Container Pull Failures

```bash
# Pre-pull containers
singularity pull docker://quay.io/biocontainers/fastp:0.23.4--h5f740d0_0

# Or use a cache directory
export NXF_SINGULARITY_CACHEDIR=/path/to/cache
```

#### Database Download Issues

```bash
# Check network connectivity
curl -I https://genome-idx.s3.amazonaws.com

# Manually download and specify custom path
--custom_kraken_db /path/to/manually/downloaded/db
```

#### SLURM Job Failures

```bash
# Check SLURM logs
cat .nextflow.log
squeue -u $USER

# Increase time limit
--max_time 480.h
```

#### Resume Not Working

```bash
# Check work directory exists
ls -la work/

# Force fresh run
nextflow run main.nf ... -resume false

# Clean and restart
nextflow clean -f
nextflow run main.nf ...
```

### Debug Mode

Enable detailed logging:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker,debug \
    -with-trace \
    -with-report \
    -with-dag
```

### Checking Logs

```bash
# Main Nextflow log
cat .nextflow.log

# Task-specific logs
cat work/*/*/.command.log
cat work/*/*/.command.err

# Find failed tasks
grep -r "ERROR" work/*/*/.command.err
```

### Resource Monitoring

```bash
# View execution report
open results/pipeline_info/execution_report_*.html

# View timeline
open results/pipeline_info/execution_timeline_*.html
```

---

## Contact and Support

- **GitHub Issues**: [github.com/gene2dis/BugBuster/issues](https://github.com/gene2dis/BugBuster/issues)
- **Email**: ffuentessantander@gmail.com
- **Documentation**: [github.com/gene2dis/BugBuster](https://github.com/gene2dis/BugBuster)

---

## Citation

If you use BugBuster in your research, please cite:

```
BugBuster: Bacterial Unraveling and metaGenomic Binning with Up-Scale Throughput, Efficient and Reproducible
Microbial Data Science Lab, Center for Bioinformatics and Integrative Biology, Universidad Andres Bello
https://github.com/gene2dis/BugBuster
```

---

*BugBuster v1.0.0 - Built with Nextflow*
