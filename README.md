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
5. Normalization of predicted genes by estimating cell number with ARGs-OAP. [`ARGs-OAP`](https://github.com/xinehc/args_oap)
6. Taxonomic profile [`Kraken2`](https://ccb.jhu.edu/software/kraken2/) or [`Sourmash`](https://sourmash.readthedocs.io/en/latest/index.html)
7. Abundance estimation [`Bracken`](https://github.com/jenniferlu717/Bracken)
8. Unification of the results with Kraken-Biom and change of format to a Phyloseq object. [`Kraken-Biom`](https://github.com/smdabdoub/kraken-biom)
9. Read traceback and taxonomic reports.
10. Genome assembly [`Megahit`](https://github.com/voutcn/megahit)
11. Contig filter [`BBmap`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
12. Taxonomic annotation of contigs using Blastn and BlobTools. [`BlobTools`](https://github.com/DRL/blobtools), [`Blast`](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
13. Functional assignation of contigs with MetaCerberus. [`MetaCerberus`](https://github.com/raw-lab/MetaCerberus)
14. ORF prediction in contigs with Prodigal. [`Prodigal`](https://github.com/hyattpd/Prodigal)
15. Prediction of resistance genes at the contig level with DeepARG. [`DeepARG`](https://github.com/gaarangoa/deeparg)
16. Contig reports, scatter plot of taxonomy at Phylum level and scatter plot of resistance genes in contigs.
17. Binning [`Metabat2`](https://bitbucket.org/berkeleylab/metabat/src/master/), [`SemiBin`](https://github.com/BigDataBiology/SemiBin) and [`COMEBin`](https://github.com/ziyewang/COMEBin)
18. Binning refinement [`MetaWrap`](https://github.com/bxlab/metaWRAP)
19. Bin quality prediction [`CheckM2`](https://github.com/chklovski/CheckM2)
20. Bin taxonomic prediction [`GTDB-TK`](https://github.com/Ecogenomics/GTDBTk)
21. Bin reports.
22. If requested functional anotation of Bins **(work in progress)** [`MetaCerberus`](https://github.com/raw-lab/MetaCerberus)
23. If requested ARG clustering **(work in progress)** [`mmseqs2`](https://github.com/soedinglab/MMseqs2)
24. Assembly modes: "coassembly", "assembly", "none"

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

All databases can be automatically download in the first use of the pipeline and their paths will be stored as symbolic links in output_path/downloaded_db folder. The descriptions of the automatic download databases are in the config/databases.config file. However, you can use your own databases by writing the absolute paths in variables prefixed with custom_ in the nextflow.config file. You can use any host bowtie2 index for filtering steps but only human host it's in the automatic download (for now 👷). 

**Note: Only the required databases for the requested tasks will be automatically download**

**Bowtie2:** must be directories with the genomes indexed with bowtie2 format
1. phiX_index = Default download from [`phage phiX174 geonome`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=genbank) (8.1 MB).
2. host_db = Default download from [`CHM13 plus Y bowtie2 index`](https://benlangmead.github.io/aws-indexes/bowtie) (4.1 GB).

**Taxonomic profiling:**
    
3A. kraken2_db = User can choice between Standard-8 (7.5 GB) and GTDB release 220 (497 GB) for automatic download. From [`kraken2 index`](https://benlangmead.github.io/aws-indexes/k2).

3B. sourmash_db = GTDB release 220 (17 GB) it's the only automatic default download in this instance (for now 👷). From [`sourmash kmers`](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs220/). 

**BBlobTools:**

4. taxdump_files = Default download from [`taxdump.tar.gz`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) (448 MB).

**BlastDB for contig taxonomic annotation:**

5. blast_db = Default download from [`ncbi nt`](https://ftp.ncbi.nlm.nih.gov/blast/db/) (434 GB).
6. deeparg_db = Default download from [`deeparg`](https://github.com/gaarangoa/deeparg) (4.8 GB).

**KARGA:**

7. karga_db = Default download from [`megares`](https://www.meglab.org/megares/download/) (9.2 MB).

**KARGVA:**

8. kargva_db = Default download from [`kargva_db`](https://github.com/DataIntellSystLab/KARGVA/tree/main) (1.5 MB).

**CheckM2**

9. checkm2_db = Default download from [`Checkm2_docs`](https://github.com/chklovski/CheckM2) (2.9 GB).

**GTDB-TK**

10. gtdbtk_db = Default download release 220 from [`gtdbtk_db`](https://ecogenomics.github.io/GTDBTk/installing/index.html) (109 GB).

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

### Feature Toggles

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--read_arg_prediction` | `false` | Read-level ARG prediction (KARGA/KARGVA) |
| `--contig_tax_and_arg` | `false` | Contig taxonomy and ARG prediction |
| `--contig_level_metacerberus` | `false` | Functional annotation with MetaCerberus |
| `--arg_bin_clustering` | `false` | ARG clustering (WIP) |

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
## Credits

gene2dis/BUGBUSTER was originally written by the Microbial Data Science Lab, Center for Bioinformatics and Integrative Biology, Universidad Andres Bello. Its development was led by Francisco A. Fuentes 

We thank the following people for their extensive assistance in the development of this pipeline:

- Francisco A. Fuentes
- Juan A. Ugalde
- Carolina Curiqueo

If you have any question of how to use the pipeline, you can contact the developer at the mail ffuentessantander@gmail.com. We will be happy to answer your questions!
