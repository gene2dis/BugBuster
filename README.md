![image](https://github.com/user-attachments/assets/a10c01f6-ef6c-40c4-a4ac-26a0c4f87564)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

## Introduction

**Bacterial Unraveling and metaGenomic Binning with Up-Scale Throughput, Efficient and Reproducible** 

**BugBuster** is a bioinformatics best-practice analysis pipeline for microbial metagenomic analysis and antimicrobial resistance gene prediction.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation uses one container per process which makes it easier to maintain and update software dependencies.

### Key Features

- **Multi-platform support**: Run locally, on HPC clusters (SLURM), or cloud (AWS, GCP, Azure)
- **Flexible execution modes**: Per-sample assembly, co-assembly, or taxonomy-only
- **Comprehensive analysis**: QC, taxonomy, assembly, binning, and ARG prediction
- **Automatic database management**: Databases download automatically on first use
- **Resource optimization**: Dynamic resource allocation with retry strategies

## Pipeline summary
![Diagrama_BugBuster drawio](https://github.com/user-attachments/assets/40d02c04-e84e-48b4-b517-c93dccd68abc)


1. Read QC, clean, and filter reads. [`FastP`](https://github.com/OpenGene/fastp)
2. Remove all samples that not have at least 10.000.000 reads after quality filter. **Can be modified by the user**
3. Remove host contamintant reads. [`Bowtie2`](https://github.com/BenLangmead/bowtie2)
4. If requested Antibiotic resistance prediction at read level using KARGA and KARGVA [`KARGA`](https://github.com/DataIntellSystLab/KARGA), [`KARGVA`](https://github.com/DataIntellSystLab/KARGVA)
5. If requested AMR gene prediction with pathogen-of-origin analysis using RGI [`RGI`](https://github.com/arpcard/rgi)
6. Normalization of predicted genes by estimating cell number with ARGs-OAP. [`ARGs-OAP`](https://github.com/xinehc/args_oap)
7. Taxonomic profile [`Kraken2`](https://ccb.jhu.edu/software/kraken2/) or [`Sourmash`](https://sourmash.readthedocs.io/en/latest/index.html)
8. Abundance estimation [`Bracken`](https://github.com/jenniferlu717/Bracken)
9. Unification of the results with Kraken-Biom and change of format to a Phyloseq object. [`Kraken-Biom`](https://github.com/smdabdoub/kraken-biom)
10. Read traceback and taxonomic reports.
11. Genome assembly [`Megahit`](https://github.com/voutcn/megahit)
12. Contig filter [`BBmap`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
13. Taxonomic annotation of contigs using Blastn and BlobTools. [`BlobTools`](https://github.com/DRL/blobtools), [`Blast`](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
14. Functional assignation of contigs with MetaCerberus. [`MetaCerberus`](https://github.com/raw-lab/MetaCerberus)
15. ORF prediction in contigs with Prodigal. [`Prodigal`](https://github.com/hyattpd/Prodigal)
16. Prediction of resistance genes at the contig level with DeepARG. [`DeepARG`](https://github.com/gaarangoa/deeparg)
17. Contig reports, scatter plot of taxonomy at Phylum level and scatter plot of resistance genes in contigs.
18. Binning with user-selectable tools (default: SemiBin; options: [`Metabat2`](https://bitbucket.org/berkeleylab/metabat/src/master/), [`SemiBin`](https://github.com/BigDataBiology/SemiBin), [`COMEBin`](https://github.com/ziyewang/COMEBin))
19. Binning refinement with [`MetaWrap`](https://github.com/bxlab/metaWRAP) (only when ≥2 binners are selected)
20. Bin quality prediction [`CheckM2`](https://github.com/chklovski/CheckM2)
21. Bin taxonomic prediction [`GTDB-TK`](https://github.com/Ecogenomics/GTDBTk)
22. Bin reports.
23. If requested functional anotation of Bins **(work in progress)** [`MetaCerberus`](https://github.com/raw-lab/MetaCerberus)
24. If requested ARG clustering **(work in progress)** [`mmseqs2`](https://github.com/soedinglab/MMseqs2)
25. Assembly modes: "coassembly", "assembly", "none"

## Quick Start

### Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=23.04.0`)
- Container runtime: [Docker](https://docs.docker.com/engine/installation/), [Singularity](https://sylabs.io/guides/), or [Podman](https://podman.io/)

### Basic Usage

```bash
# Clone the repository
git clone https://github.com/gene2dis/BugBuster.git
cd BugBuster

# Run with Docker
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker

# Run with Singularity (HPC)
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile singularity
```

### Available Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers |
| `podman` | Run with Podman containers |
| `slurm` | Run on SLURM HPC cluster |
| `aws` | Run on AWS Batch |
| `gcp` | Run on Google Cloud |
| `azure` | Run on Azure Batch |
| `test` | Run with minimal test data |

Combine profiles as needed: `-profile slurm,singularity` or `-profile aws,docker`

### Cloud Execution Examples

```bash
# AWS Batch
nextflow run main.nf \
    --input s3://bucket/samplesheet.csv \
    --output s3://bucket/results \
    -profile aws,docker \
    -work-dir s3://bucket/work

# Google Cloud
nextflow run main.nf \
    --input gs://bucket/samplesheet.csv \
    --output gs://bucket/results \
    -profile gcp,docker \
    --gcp_project my-project \
    -work-dir gs://bucket/work
```

See [docs/deployment.md](docs/deployment.md) for detailed deployment instructions.

## Databases

All databases are automatically downloaded on first use and stored in `<output_dir>/../databases/` by default (configurable via `--databases_dir`). Symbolic links are created in `<output_dir>/downloaded_db/` for reference. Only databases required for your selected pipeline options will be downloaded.

You can use custom databases by specifying paths with `--custom_*` parameters (see table below).

### Automatic Download Databases

| Database | Size | Used For | Trigger Parameter | Custom Path Parameter |
|----------|------|----------|-------------------|----------------------|
| **phiX174 Index** | 8.1 MB | PhiX contamination removal | `quality_control=true` | `--custom_phiX_index` |
| **Human Host Index** | 4.1 GB | Host read removal | `quality_control=true` | `--custom_bowtie_host_index` |
| **Kraken2 Standard-8** | 7.5 GB | Taxonomic profiling | `taxonomic_profiler='kraken2'` | `--custom_kraken_db` |
| **Kraken2 GTDB r220** | 497 GB | Taxonomic profiling | `kraken2_db='gtdb_220'` | `--custom_kraken_db` |
| **Sourmash GTDB r220** | 17 GB | Taxonomic profiling | `taxonomic_profiler='sourmash'` | `--custom_sourmash_db` |
| **NCBI Taxdump** | 448 MB | BlobTools taxonomy | `contig_tax_and_arg=true` | `--custom_taxdump_files` |
| **NCBI NT** | 434 GB | Contig BLAST | `contig_tax_and_arg=true` | `--custom_blast_db` |
| **DeepARG** | 4.8 GB | Contig ARG prediction | `contig_tax_and_arg=true` | `--custom_deeparg_db` |
| **KARGA (MEGARes)** | 9.2 MB | Read ARG prediction | `read_arg_prediction=true` | `--custom_karga_db` |
| **KARGVA** | 1.5 MB | Read ARG variant prediction | `read_arg_prediction=true` | `--custom_kargva_db` |
| **CARD (RGI)** | 500 MB - 50 GB | AMR gene prediction with pathogen-of-origin | `rgi_prediction=true` | `--custom_rgi_card_db`, `--custom_rgi_wildcard` |
| **CheckM2** | 2.9 GB | Bin quality assessment | `include_binning=true` | `--custom_checkm2_db` |
| **GTDB-TK r220** | 109 GB | Bin taxonomic classification | `include_binning=true` | `--custom_gtdbtk_db` |

### Database Sources

- **phiX_index**: [`phage phiX174 genome`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=genbank)
- **host_db**: [`CHM13 plus Y bowtie2 index`](https://benlangmead.github.io/aws-indexes/bowtie)
- **kraken2_db**: [`kraken2 index`](https://benlangmead.github.io/aws-indexes/k2)
- **sourmash_db**: [`sourmash kmers`](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs220/)
- **taxdump_files**: [`taxdump.tar.gz`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
- **blast_db**: [`ncbi nt`](https://ftp.ncbi.nlm.nih.gov/blast/db/)
- **deeparg_db**: [`deeparg`](https://github.com/gaarangoa/deeparg)
- **karga_db**: [`megares`](https://www.meglab.org/megares/download/)
- **kargva_db**: [`kargva_db`](https://github.com/DataIntellSystLab/KARGVA/tree/main)
- **card_db**: [`CARD`](https://card.mcmaster.ca/) - Comprehensive Antibiotic Resistance Database
- **checkm2_db**: [`Checkm2_docs`](https://github.com/chklovski/CheckM2)
- **gtdbtk_db**: [`gtdbtk_db`](https://ecogenomics.github.io/GTDBTk/installing/index.html)

## Samplesheet Format

Create a CSV file with your sample information:

```csv
sample,r1,r2,s
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,/path/to/sample2_singletons.fastq.gz
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `r1` | Yes | Path to forward reads (R1) |
| `r2` | Yes | Path to reverse reads (R2) |
| `s` | No | Path to singleton reads (optional) |

## Pipeline Parameters

### Core Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | *required* | Path to samplesheet CSV |
| `--output` | *required* | Output directory |
| `--quality_control` | `true` | Enable QC and host filtering |
| `--assembly_mode` | `assembly` | `assembly`, `coassembly`, or `none` |
| `--taxonomic_profiler` | `sourmash` | `kraken2`, `sourmash`, or `none` |
| `--include_binning` | `false` | Enable binning and refinement |
| `--binners` | `semibin` | Binners to run: `comebin`, `semibin`, `metabat2` (comma-separated; ≥2 enables MetaWRAP) |
| `--min_read_sample` | `0` | Minimum reads required after QC |

### Feature Toggles

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--read_arg_prediction` | `false` | Read-level ARG prediction (KARGA/KARGVA) |
| `--rgi_prediction` | `false` | AMR gene prediction with pathogen-of-origin (RGI/CARD) |
| `--contig_tax_and_arg` | `false` | Contig taxonomy and ARG prediction |
| `--contig_level_metacerberus` | `false` | Functional annotation with MetaCerberus |
| `--arg_bin_clustering` | `false` | ARG clustering (WIP) |

### Database Selection

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `--kraken2_db` | `standard-8` | `standard-8`, `gtdb_220` | Kraken2 database version |
| `--sourmash_db` | `gtdb_220_k31` | `gtdb_220_k31` | Sourmash database version |
| `--databases_dir` | `<output>/../databases` | Path | Database storage location |

### Quality Control Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastp_qualified_quality_phred` | `20` | Phred quality threshold |
| `--fastp_unqualified_percent_limit` | `10` | Max % unqualified bases |
| `--fastp_n_base_limit` | `5` | Max N bases per read |

### Taxonomy Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--kraken_confidence` | `0.1` | Kraken2 confidence threshold (0-1) |
| `--bracken_tax_level` | `S` | Bracken level: D, P, C, O, F, G, S |
| `--sourmash_tax_rank` | `species` | Sourmash rank: genus, species, strain |
| `--taxonomy_plot_levels` | `Phylum,Family,Genus,Species` | Taxonomic levels to plot |
| `--taxonomy_top_n_taxa` | `10` | Number of top taxa in plots |
| `--create_phyloseq_rds` | `false` | Generate R phyloseq RDS files |

### Assembly & Binning Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bbmap_lenght` | `1000` | Minimum contig length after filtering |
| `--binners` | `semibin` | Binners to run (comma-separated). `≥2` selected → MetaWRAP refinement enabled |
| `--metabat_minContig` | `2500` | Minimum contig length for MetaBAT2 |
| `--metawrap_completeness` | `50` | Minimum bin completeness (%) — used only when ≥2 binners |
| `--metawrap_contamination` | `10` | Maximum bin contamination (%) — used only when ≥2 binners |
| `--semibin_env_model` | `human_gut` | SemiBin environment model |

**SemiBin Models**: `human_gut`, `dog_gut`, `cat_gut`, `mouse_gut`, `pig_gut`, `human_oral`, `chicken_caecum`, `ocean`, `soil`, `wastewater`, `built_environment`, `global`

**Binner selection examples:**
```bash
# Default: single fast binner, no MetaWRAP
--include_binning true --binners semibin

# Two binners + MetaWRAP refinement
--include_binning true --binners semibin,metabat2

# All three binners + MetaWRAP
--include_binning true --binners semibin,metabat2,comebin
```

### ARG Prediction Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--deeparg_min_prob` | `0.8` | Minimum probability threshold |
| `--deeparg_arg_alignment_identity` | `50` | Minimum alignment identity (%) |
| `--deeparg_arg_alignment_evalue` | `1e-10` | Maximum E-value |
| `--deeparg_model_version` | `v2` | DeepARG model version |
| `--rgi_card_version` | `latest` | CARD database version |
| `--rgi_include_wildcard` | `true` | Include WildCARD variants |
| `--rgi_aligner` | `kma` | RGI aligner: kma, bowtie2, bwa |
| `--rgi_kmer_size` | `61` | K-mer size for pathogen prediction |
| `--rgi_min_kmer_coverage` | `10` | Minimum k-mer coverage |

### Resource Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `240.h` | Maximum time per process |

### Help

```bash
nextflow run main.nf --help
```

## Output Structure

The pipeline generates organized output with numbered prefixes:

```
results/
├── pipeline_info/              # Execution reports and logs
├── 01_quality_control/         # QC results (if quality_control=true)
│   ├── fastp/                  # Per-sample FastP reports
│   └── summary/                # Aggregated statistics
├── 02_taxonomy/                # Taxonomic profiling (if taxonomic_profiler != 'none')
│   ├── kraken2/ or sourmash/   # Per-sample results
│   ├── tables/                 # Abundance tables
│   ├── phyloseq/               # Phyloseq objects (*.h5, *.RDS)
│   └── figures/                # Taxonomy plots
├── 03_assembly/                # Genome assembly (if assembly_mode != 'none')
│   ├── per_sample/             # Per-sample assemblies
│   └── coassembly/             # Co-assembly results
├── 04_binning/                 # Metagenomic binning (if include_binning=true)
│   ├── per_sample/ or coassembly/
│   │   ├── raw_bins/           # MetaBAT2, SemiBin, COMEBin
│   │   ├── refined_bins/       # MetaWRAP refined bins
│   │   ├── quality/            # CheckM2 reports
│   │   └── taxonomy/           # GTDB-TK classifications
│   ├── quality/summary/        # Aggregated quality reports
│   └── taxonomy/summary/       # Aggregated taxonomy reports
├── 05_arg_prediction/          # ARG predictions
│   ├── read_level/             # Read-level predictions
│   │   ├── karga/              # KARGA results (if read_arg_prediction=true)
│   │   ├── kargva/             # KARGVA results (if read_arg_prediction=true)
│   │   ├── rgi/                # RGI per-sample results (if rgi_prediction=true)
│   │   ├── rgi_kmer/           # RGI pathogen-of-origin (if rgi_prediction=true)
│   │   └── rgi_summary/        # RGI aggregated reports (if rgi_prediction=true)
│   ├── contig_level/           # DeepARG (if contig_tax_and_arg=true)
│   └── bin_level/              # Bin ARG clustering (if arg_bin_clustering=true)
├── 06_contig_taxonomy/         # BlobTools plots (if contig_tax_and_arg=true)
└── 07_functional_annotation/   # MetaCerberus (if contig_level_metacerberus=true)
```

**Database Storage**: Databases are stored separately at `<output_dir>/../databases/` by default (configurable via `--databases_dir`).

For detailed output descriptions, see [`docs/manual.md`](docs/manual.md#8-output-structure).

## Credits

gene2dis/BUGBUSTER was originally written by the Microbial Data Science Lab, Center for Bioinformatics and Integrative Biology, Universidad Andres Bello. Its development was led by Francisco A. Fuentes 

We thank the following people for their extensive assistance in the development of this pipeline:

- Francisco A. Fuentes
- Juan A. Ugalde
- Carolina Curiqueo

If you have any question of how to use the pipeline, you can contact the developer at the mail ffuentessantander@gmail.com. We will be happy to answer your questions!
