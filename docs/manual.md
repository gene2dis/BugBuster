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
| **CARD (RGI)** | 500 MB - 50 GB | AMR gene prediction with pathogen-of-origin | `rgi_prediction=true` |
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
| `--binners` | `semibin` | `comebin`, `semibin`, `metabat2` | Comma-separated binners to run; ≥2 enables MetaWRAP |
| `--read_arg_prediction` | `false` | `true`, `false` | Read-level ARG prediction (KARGA/KARGVA) |
| `--rgi_prediction` | `false` | `true`, `false` | AMR prediction with pathogen-of-origin (RGI/CARD) |
| `--contig_tax_and_arg` | `false` | `true`, `false` | Contig taxonomy and ARG prediction |
| `--metacerberus_levels` | `none` | `reads`, `contigs`, `bins`, `none` | Comma-separated levels to run MetaCerberus functional annotation |
| `--metacerberus_hmm` | `KOFam_all,COG,VOG,PHROG,CAZy` | Any valid HMM database names | Comma-separated HMM databases for MetaCerberus |
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
| `--custom_rgi_card_db` | Path to pre-prepared CARD database directory |
| `--custom_rgi_wildcard` | Path to WildCARD directory (use with `--custom_rgi_card_db`) |

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
| `--bbmap_lenght` | `1000` | Minimum contig length after BBMap filtering |

### 6.8 Binning Options

#### Binner Selection

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `--binners` | `semibin` | `comebin`, `semibin`, `metabat2` | Comma-separated list of binners to run. At least one required. MetaWRAP is automatically used when ≥2 binners are selected. |

**Selection logic:**
- **1 binner selected** → runs that binner only, MetaWRAP skipped, binner output used directly
- **≥2 binners selected** → runs selected binners + MetaWRAP bin refinement

#### Basic Binning Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metabat_minContig` | `2500` | Minimum contig length for MetaBAT2 |
| `--metawrap_completeness` | `50` | Minimum bin completeness (%) for MetaWRAP (used when ≥2 binners) |
| `--metawrap_contamination` | `10` | Maximum bin contamination (%) for MetaWRAP (used when ≥2 binners) |
| `--semibin_env_model` | `human_gut` | SemiBin environment model |

**SemiBin Environment Models:**
- `human_gut`, `dog_gut`, `cat_gut`, `mouse_gut`, `pig_gut`
- `human_oral`, `chicken_caecum`
- `ocean`, `soil`, `wastewater`, `built_environment`
- `global` (for mixed/unknown environments)

#### Advanced MetaBAT2 Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metabat_maxP` | `95` | Maximum percentage of good contigs |
| `--metabat_minS` | `60` | Minimum score for binning |
| `--metabat_maxEdges` | `200` | Maximum edges in the graph |
| `--metabat_pTNF` | `0` | TNF probability threshold |
| `--metabat_minCV` | `1` | Minimum coefficient of variation |
| `--metabat_minCVSum` | `1` | Minimum sum of coefficient of variation |
| `--metabat_minClsSize` | `200000` | Minimum cluster size (bp) |

### 6.9 DeepARG Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--deeparg_min_prob` | `0.8` | Minimum probability threshold (0-1) |
| `--deeparg_arg_alignment_identity` | `50` | Minimum alignment identity (%) |
| `--deeparg_arg_alignment_evalue` | `1e-10` | Maximum E-value |
| `--deeparg_arg_alignment_overlap` | `0.8` | Minimum alignment overlap (0-1) |
| `--deeparg_arg_num_alignments_per_entry` | `1000` | Number of alignments per entry |
| `--deeparg_model_version` | `v2` | DeepARG model version (`v1` or `v2`) |

### 6.10 RGI Options

Parameters for RGI AMR gene prediction with pathogen-of-origin analysis:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--rgi_card_version` | `latest` | CARD database version (`latest` or specific version) |
| `--rgi_include_wildcard` | `true` | Include WildCARD variants for extended allelic diversity |
| `--rgi_aligner` | `kma` | Read aligner: `kma`, `bowtie2`, or `bwa` |
| `--rgi_kmer_size` | `61` | K-mer size for pathogen-of-origin prediction |
| `--rgi_min_kmer_coverage` | `10` | Minimum k-mer coverage threshold |

**Note**: For detailed manual CARD database preparation instructions, see [`docs/RGI_IMPLEMENTATION_PLAN.md`](RGI_IMPLEMENTATION_PLAN.md).

### 6.11 Bowtie2 Alignment Options

Advanced parameters for Bowtie2 read alignment during host filtering:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bowtie_ma` | `2` | Match bonus |
| `--bowtie_mp` | `6,2` | Mismatch penalty (max, min) |
| `--bowtie_score_min` | `G,15,6` | Minimum alignment score function |
| `--bowtie_k` | `1` | Number of alignments to report |
| `--bowtie_N` | `1` | Number of mismatches allowed in seed |
| `--bowtie_L` | `20` | Seed length |
| `--bowtie_R` | `2` | Number of re-seeding attempts |
| `--bowtie_i` | `S,1,0.75` | Interval function for seeding |

### 6.12 MMseqs2 Clustering Options

Parameters for ARG clustering (used when `--arg_bin_clustering=true`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mmseqs_start_sens` | `2` | Starting sensitivity |
| `--mmseqs_s` | `7` | Sensitivity level |
| `--mmseqs_sens_steps` | `3` | Number of sensitivity steps |
| `--mmseqs_min_seq_id` | `0.8` | Minimum sequence identity (0-1) |
| `--mmseqs_c` | `0.7` | Coverage threshold (0-1) |
| `--mmseqs_cov_mode` | `2` | Coverage mode |
| `--mmseqs_e` | `1e-20` | E-value threshold |
| `--mmseqs_format_mode` | `4` | Output format mode |
| `--mmseqs_alignment_mode` | `3` | Alignment mode |
| `--mmseqs_max_seqs` | `10000` | Maximum number of sequences |
| `--mmseqs_format_output` | `empty,query,target,evalue,pident,qcov,tcov,tseq` | Output format fields |

### 6.12 MetaCerberus Options

MetaCerberus is **disabled by default** (`--metacerberus_levels none`). Activate it by specifying one or more levels.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metacerberus_levels` | `none` | Comma-separated levels to run: `reads`, `contigs`, `bins`, or `none` |
| `--metacerberus_hmm` | `KOFam_all,KOFam_prokaryote,COG,VOG,PHROG,CAZy,MetHMMDB,Pfam` | Comma-separated HMM databases to use |
| `--metacerberus_db` | `<databases_dir>/metacerberus` | Path to directory where HMM databases are stored |
| `--metacerberus_minscore` | `25` | Minimum HMM score |
| `--metacerberus_evalue` | `1e-09` | Maximum E-value |

**Available HMM Databases:** `KOFam_all`, `KOFam_eukaryote`, `KOFam_prokaryote`, `COG`, `VOG`, `PHROG`, `CAZy`, `MetHMMDB`, `Pfam`

### 6.13 Taxonomy Visualization Options

Control taxonomic output visualization and formatting:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--taxonomy_plot_levels` | `Phylum,Family,Genus,Species` | Comma-separated taxonomic levels to plot |
| `--taxonomy_top_n_taxa` | `10` | Number of top taxa to display in plots |
| `--create_phyloseq_rds` | `false` | Generate R phyloseq RDS files for downstream analysis |

**Taxonomic Levels:** Domain (D), Phylum (P), Class (C), Order (O), Family (F), Genus (G), Species (S)

### 6.14 Database Storage Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--databases_dir` | `<output>/../databases` | Directory for storing downloaded databases (separate from results) |

By default, databases are stored in a `databases/` directory at the same level as your output directory. This allows database reuse across multiple pipeline runs. Symbolic links are created in `<output>/downloaded_db/` for reference.

### 6.15 Resource Limit Options

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

### Main Output Directory

```
results/
├── pipeline_info/                              # Execution reports
│   ├── execution_report_*.html                 # Resource usage report
│   ├── execution_timeline_*.html               # Timeline visualization
│   ├── execution_trace_*.txt                   # Task trace log
│   └── pipeline_dag_*.svg                      # Pipeline DAG
├── 01_quality_control/                         # Quality control (if quality_control=true)
│   ├── fastp/                                  # FastP reports per sample
│   │   └── {sample}/                           # Per-sample QC results
│   │       ├── {sample}.fastp.html             # HTML report
│   │       ├── {sample}.fastp.json             # JSON report
│   │       └── {sample}.fastp.log              # Log file
│   └── summary/                                # Aggregated QC statistics
│       ├── Reads_report.csv                    # Read count summary
│       └── *.png                               # QC plots
├── 02_taxonomy/                                # Taxonomic profiling (if taxonomic_profiler != 'none')
│   ├── kraken2/                                # Kraken2 results (if taxonomic_profiler='kraken2')
│   │   └── {sample}_*.report.txt               # Per-sample Kraken2 reports
│   ├── bracken/                                # Bracken abundance estimates
│   │   └── {sample}_*.bracken                  # Per-sample Bracken results
│   ├── sourmash/                               # Sourmash results (if taxonomic_profiler='sourmash')
│   │   └── *.with-lineages.csv                 # Per-sample Sourmash results
│   ├── tables/                                 # Taxonomy tables
│   │   └── *.tsv                               # Abundance tables
│   ├── phyloseq/                               # Phyloseq objects
│   │   ├── *.h5                                # HDF5 format
│   │   └── *.RDS                               # R phyloseq object (if create_phyloseq_rds=true)
│   └── figures/                                # Taxonomy plots
│       └── *.png                               # Visualization plots
├── 03_assembly/                                # Genome assembly (if assembly_mode != 'none')
│   ├── per_sample/                             # Per-sample assemblies (if assembly_mode='assembly')
│   │   └── {sample}/
│   │       ├── {sample}_filtered_contigs.fa    # Filtered contigs (≥ bbmap_lenght bp)
│   │       ├── {sample}_contig.stats           # Assembly statistics
│   │       └── {sample}_contigs.fa             # Raw MEGAHIT contigs
│   └── coassembly/                             # Co-assembly (if assembly_mode='coassembly')
│       ├── coassembly_filtered_contigs.fa      # Filtered contigs
│       ├── coassembly_contig.stats             # Assembly statistics
│       └── coassembly_contigs.fa               # Raw MEGAHIT contigs
├── 04_binning/                                 # Metagenomic binning (if include_binning=true)
│   ├── per_sample/                             # Per-sample binning (if assembly_mode='assembly')
│   │   └── {sample}/
│   │       ├── raw_bins/                       # Raw bins from 3 binners
│   │       │   ├── metabat2/                   # MetaBAT2 bins
│   │       │   ├── semibin/                    # SemiBin bins
│   │       │   └── comebin/                    # COMEBin bins
│   │       ├── refined_bins/                   # MetaWRAP refined bins
│   │       ├── quality/                        # CheckM2 quality reports
│   │       │   └── checkm2/
│   │       └── taxonomy/                       # GTDB-TK taxonomy
│   │           └── gtdbtk/
│   ├── coassembly/                             # Co-assembly binning (if assembly_mode='coassembly')
│   │   ├── raw_bins/                           # Raw bins from 3 binners
│   │   ├── refined_bins/                       # MetaWRAP refined bins
│   │   ├── quality/                            # CheckM2 quality reports
│   │   ├── taxonomy/                           # GTDB-TK taxonomy
│   │   ├── coverage/                           # Bin coverage information
│   │   └── summary/                            # Bin summary statistics
│   ├── quality/                                # Aggregated quality reports
│   │   └── summary/
│   │       ├── *.csv                           # Quality summary tables
│   │       └── *.png                           # Quality plots
│   └── taxonomy/                               # Aggregated taxonomy reports
│       └── summary/
│           ├── *.csv                           # Taxonomy summary tables
│           └── *.png                           # Taxonomy plots
├── 05_arg_prediction/                          # ARG predictions
│   ├── read_level/                             # Read-level ARG
│   │   ├── karga/                              # KARGA results (if read_arg_prediction=true)
│   │   │   └── {sample}_KARGA_mappedGenes.csv
│   │   ├── kargva/                             # KARGVA results (if read_arg_prediction=true)
│   │   │   └── {sample}_KARGVA_mappedGenes.csv
│   │   ├── args_oap/                           # ARGs-OAP normalization (if read_arg_prediction=true)
│   │   │   └── {sample}_args_oap_s1_out/
│   │   ├── rgi/                                # RGI per-sample results (if rgi_prediction=true)
│   │   │   └── {sample}/
│   │   │       ├── *_allele_mapping_data.txt
│   │   │       ├── *_gene_mapping_data.txt
│   │   │       └── *_sorted.length_100.bam
│   │   ├── rgi_kmer/                           # RGI pathogen-of-origin (if rgi_prediction=true)
│   │   │   └── {sample}/
│   │   │       ├── *_61mer_analysis.json
│   │   │       └── *_61mer_analysis.txt
│   │   ├── rgi_summary/                        # RGI aggregated reports (if rgi_prediction=true)
│   │   │   ├── RGI_summary_report.csv
│   │   │   ├── RGI_amr_gene_family_distribution.png
│   │   │   ├── RGI_drug_class_profile.png
│   │   │   ├── RGI_resistance_mechanisms.png
│   │   │   └── RGI_sample_amr_counts.png
│   │   └── summary/                            # Normalized ARG summary (if read_arg_prediction=true)
│   │       └── *.csv
│   ├── contig_level/                           # Contig-level ARG (if contig_tax_and_arg=true)
│   │   ├── prodigal/                           # ORF predictions
│   │   │   └── {sample}/
│   │   ├── deeparg/                            # DeepARG predictions (not published by default)
│   │   ├── summary/                            # ARG summary reports
│   │   │   └── Contig_tax_and_arg_prediction.tsv
│   │   └── figures/                            # ARG visualization
│   │       └── *.png
│   └── bin_level/                              # Bin-level ARG (if arg_bin_clustering=true)
│       ├── proteins/                           # Prodigal ORF predictions
│       │   └── {sample}/
│       ├── deeparg/                            # DeepARG predictions per bin
│       │   └── {sample}/
│       └── clustering/                         # MMseqs2 clustering results
│           └── *_cluster.tsv
├── 06_contig_taxonomy/                         # Contig taxonomy (if contig_tax_and_arg=true)
│   └── figures/                                # BlobTools plots
│       └── *.png
└── 07_functional_annotation/                   # Functional annotation (if metacerberus_levels!=none)
    ├── reads/                                  # Read-level annotation
    │   └── {sample}/
    │       └── {sample}_annotation_results/
    ├── contigs/                                # Contig-level annotation
    │   └── {sample}/
    │       └── {sample}_annotation_results/
    └── bins/                                   # Bin-level annotation
        └── {sample}_{bin}_annotation_results/
```

### Database Storage Directory

By default, databases are stored separately from results at `<output_dir>/../databases/`:

```
databases/                            # Database storage (configurable via --databases_dir)
├── phiX_index/                       # PhiX174 Bowtie2 index
├── host_index/                       # Host genome Bowtie2 index
├── kraken2_standard8/                # Kraken2 Standard-8 database
├── kraken2_gtdb220/                  # Kraken2 GTDB r220 database (if selected)
├── sourmash_gtdb220/                 # Sourmash GTDB r220 database
├── taxdump/                          # NCBI taxonomy dump
├── blast_nt/                         # NCBI NT BLAST database
├── deeparg/                          # DeepARG database
├── karga/                            # KARGA (MEGARes) database
├── kargva/                           # KARGVA database
├── rgi/                              # CARD database for RGI
├── checkm2/                          # CheckM2 database
└── gtdbtk_r220/                      # GTDB-TK release 220 database
```

### Conditional Outputs

The following outputs are only generated when specific parameters are enabled:

| Output Directory | Required Parameter | Description |
|------------------|-------------------|-------------|
| `01_quality_control/` | `quality_control=true` | Quality control and filtering results |
| `02_taxonomy/` | `taxonomic_profiler != 'none'` | Taxonomic profiling results |
| `02_taxonomy/phyloseq/*.RDS` | `create_phyloseq_rds=true` | R phyloseq object for downstream analysis |
| `03_assembly/` | `assembly_mode != 'none'` | Assembly results |
| `03_assembly/per_sample/` | `assembly_mode='assembly'` | Per-sample assemblies |
| `03_assembly/coassembly/` | `assembly_mode='coassembly'` | Co-assembly results |
| `04_binning/` | `include_binning=true` | Binning and bin refinement results |
| `05_arg_prediction/read_level/karga/` | `read_arg_prediction=true` | KARGA/KARGVA ARG predictions |
| `05_arg_prediction/read_level/rgi/` | `rgi_prediction=true` | RGI AMR predictions with pathogen-of-origin |
| `05_arg_prediction/contig_level/` | `contig_tax_and_arg=true` | Contig-level ARG predictions |
| `05_arg_prediction/bin_level/` | `arg_bin_clustering=true` | Bin-level ARG clustering |
| `06_contig_taxonomy/` | `contig_tax_and_arg=true` | Contig taxonomic annotation (BlobTools) |
| `07_functional_annotation/` | `metacerberus_levels!=none` | Functional annotation results (reads/contigs/bins) |

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

### Example 3: Full Pipeline with Binning (Single Binner - Fast)

Complete analysis using the default single binner (no MetaWRAP overhead):

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --quality_control true \
    --assembly_mode assembly \
    --taxonomic_profiler sourmash \
    --include_binning true \
    --binners semibin \
    --semibin_env_model human_gut \
    -profile singularity
```

### Example 3b: Full Pipeline with Binning + MetaWRAP Refinement

Run multiple binners and combine results with MetaWRAP for higher quality MAGs:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --quality_control true \
    --assembly_mode assembly \
    --taxonomic_profiler sourmash \
    --include_binning true \
    --binners semibin,metabat2 \
    --semibin_env_model human_gut \
    --metawrap_completeness 50 \
    --metawrap_contamination 10 \
    -profile singularity
```

**YAML equivalent (`params.yaml`):**
```yaml
input: "samplesheet.csv"
output: "./results"
quality_control: true
assembly_mode: "assembly"
taxonomic_profiler: "sourmash"
include_binning: true
binners: "semibin,metabat2"
semibin_env_model: "human_gut"
metawrap_completeness: 50
metawrap_contamination: 10
```
```bash
nextflow run main.nf -params-file params.yaml -profile singularity
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
    --rgi_prediction true \
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

### Example 11: Custom Taxonomy Visualization

Enhanced taxonomy profiling with custom visualization:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --taxonomic_profiler kraken2 \
    --kraken2_db gtdb_220 \
    --taxonomy_plot_levels "Phylum,Class,Order,Family,Genus" \
    --taxonomy_top_n_taxa 20 \
    --create_phyloseq_rds true \
    --assembly_mode none \
    -profile docker
```

### Example 12: Shared Database Storage

Using pre-downloaded databases in shared storage:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --databases_dir /shared/databases/bugbuster \
    --custom_kraken_db /shared/databases/bugbuster/kraken2_gtdb220 \
    --custom_checkm2_db /shared/databases/bugbuster/checkm2/uniref100.KO.1.dmnd \
    --custom_gtdbtk_db /shared/databases/bugbuster/gtdbtk_r220 \
    --include_binning true \
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
