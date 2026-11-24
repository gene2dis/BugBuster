![image](https://github.com/user-attachments/assets/a10c01f6-ef6c-40c4-a4ac-26a0c4f87564)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

[![BugBuster overview: reads go through quality control, then read-level taxonomy, ARG and functional profiling, assembly with binning into MAGs, and contig-level ARG, taxonomy, mobile-element and functional branches; every stage feeds the run summary. Green boxes run by default, dashed boxes are optional.](docs/images/bugbuster_overview.svg)](docs/images/bugbuster_overview.png)

Green boxes run by default; dashed boxes are switched on by a parameter. Every stage and its switch: [docs/WORKFLOW_DIAGRAM.md](docs/WORKFLOW_DIAGRAM.md).

## Introduction

**Bacterial Unraveling and metaGenomic Binning with Up-Scale Throughput, Efficient and Reproducible**

**BugBuster** is a [Nextflow](https://www.nextflow.io) pipeline for shotgun metagenomes: read quality control and host removal, taxonomic and functional profiles of the reads, assembly and binning into genomes, and antimicrobial resistance genes (ARGs) at the read, contig and genome level. Every step runs in a container (Docker, Singularity/Apptainer or Podman), on a workstation, a SLURM cluster or AWS, Google Cloud and Azure, and the run ends with one summary of every stage that ran.

### Key Features

- **Runs anywhere Nextflow does**: locally, on SLURM, or on AWS Batch, Google Cloud Batch and Azure Batch
- **Per-sample assembly, co-assembly (also per metadata group), or no assembly**
- **ARGs three ways**: from the reads (RGI with pathogen of origin, K-MARVEL, ARGs-OAP), on contigs with their taxonomy and mobile-element class, and per genome
- **Databases managed for you**: downloaded on first use, stored once per `--databases_dir` and reused by later runs
- **One results summary**: tables, an Excel workbook, figures, an offline HTML report and a PDF report

## Pipeline at a glance

Stages in the order they run; the default column says which run without any option.

| Stage | Tools | Switch | Default |
|-------|-------|--------|---------|
| Read quality control | [fastp](https://github.com/OpenGene/fastp) | `--quality_control` | on |
| Host and phiX removal | [bwa-mem3](https://github.com/fg-labs/bwa-mem3), or [Bowtie2](https://github.com/BenLangmead/bowtie2) with `--read_aligner bowtie2` | `--quality_control` (with fastp) | on, bwa-mem3 |
| Sequencing coverage | [Nonpareil](https://github.com/lmrodriguezr/nonpareil) | `--nonpareil` | on |
| Taxonomic profile | [sylph](https://github.com/bluenote-1577/sylph) + [sylph-tax](https://github.com/bluenote-1577/sylph-tax), [Kraken2](https://ccb.jhu.edu/software/kraken2/) + [Bracken](https://github.com/jenniferlu717/Bracken), or [sourmash](https://sourmash.readthedocs.io/en/latest/index.html) | `--taxonomic_profiler` | sylph |
| Read-level ARGs | [RGI](https://github.com/arpcard/rgi) (with pathogen of origin), [K-MARVEL](https://bitbucket.org/amr-avengers/k-marvel), [ARGs-OAP](https://github.com/xinehc/args_oap) | `--rgi_prediction`, `--kmarvel_prediction`, `--args_oap_prediction` | off |
| Read-level functions | [Woltka](https://github.com/qiyunzhu/woltka), [SUPER-FOCUS](https://github.com/metageni/SUPER-FOCUS) or [HUMAnN 4](https://github.com/biobakery/humann) | `--read_level_functional` | off |
| Assembly | [MEGAHIT](https://github.com/voutcn/megahit), contig length filter with [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) | `--assembly_mode` | on, per sample |
| Binning into MAGs | [SemiBin2](https://github.com/BigDataBiology/SemiBin), [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/), [COMEBin](https://github.com/ziyewang/COMEBin); [Binette](https://github.com/genotoul-bioinfo/Binette) refinement with two or more; [CheckM2](https://github.com/chklovski/CheckM2), [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | `--include_binning`, `--binners` | off |
| MAG checks | [GUNC](https://github.com/grp-bork/gunc) chimerism; multi-sample binning; dereplication with [galah](https://github.com/wwood/galah) + [CoverM](https://github.com/wwood/CoverM) | `--gunc_prediction`, `--multisample_binning`, `--dereplicate_mags` | off |
| Contig taxonomy and ARGs | [Pyrodigal](https://github.com/althonos/pyrodigal), [Metabuli](https://github.com/steineggerlab/Metabuli), [BlobTools](https://github.com/DRL/blobtools), [DeepARG](https://github.com/gaarangoa/deeparg) | `--contig_tax_and_arg` | off |
| Plasmids and viruses | [geNomad](https://github.com/apcamargo/genomad), [CheckV](https://bitbucket.org/berkeleylab/checkv) | `--genomad_prediction`, `--checkv_prediction` | off |
| Contig functions | [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper), [run_dbcan](https://github.com/bcb-unl/run_dbcan), [featureCounts](https://subread.sourceforge.net), [MicrobeCensus](https://github.com/snayfach/MicrobeCensus) | `--contig_level_functional` | off |
| MAG annotation | [Bakta](https://github.com/oschwengers/bakta) | `--mag_level_functional` | off |
| ARG clusters across MAGs | [MMseqs2](https://github.com/soedinglab/MMseqs2) | `--arg_bin_clustering` | off |
| Results summary | tables, workbook, figures, HTML and PDF report | `--results_summary` | on |

What each stage needs and costs: [manual, section 5](docs/manual.md#5-choosing-what-runs); every stage with the parameter that gates it: [docs/WORKFLOW_DIAGRAM.md](docs/WORKFLOW_DIAGRAM.md).

Terms used throughout: **singleton reads** are reads whose mate was lost (for example in an earlier trimming step), given in the samplesheet's optional `s` column; a **MAG** (metagenome-assembled genome) is a genome bin recovered from an assembly; **CPGE** (copies per genome equivalent) is a function's read abundance divided by the sample's genome equivalents, the number of genomes its reads amount to (estimated by MicrobeCensus).

## Quick Start

### Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) 25.10.0 or later, with Java 17 to 25 (17 to 26 with Nextflow 26.04)
- A container engine: [Docker](https://docs.docker.com/engine/installation/), [Apptainer](https://apptainer.org/docs/admin/main/installation.html) or Singularity 3.6+, or [Podman](https://podman.io/)
- For the default aligner, bwa-mem3: an x86-64 CPU with AVX2 (Haswell 2013 or newer), and ~75 GB of RAM once, to build its human index. With less RAM or an older CPU, use `--read_aligner bowtie2`. Other stages' needs: [hardware requirements](docs/manual.md#hardware-requirements)

### Install and check the installation

```bash
git clone https://github.com/gene2dis/BugBuster.git
cd BugBuster

# Test run: two small samples from GitHub, QC only, no database download
nextflow run main.nf -profile test,docker    # or test,singularity / test,apptainer
```

The test run takes a few minutes (about two once the container images are pulled), writes its results to `test_results/` and ends with `Pipeline completed successfully`. What it runs: [manual, section 3](docs/manual.md#check-the-installation-test-profile).

### First real run

Write a samplesheet (template: `assets/samplesheet_example.csv`; paired-end reads only; rules for names and paths in [the manual](docs/manual.md#samplesheet-format)):

```csv
sample,r1,r2,s
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,/path/to/sample2_singletons.fastq.gz
```

```bash
nextflow run main.nf --input samplesheet.csv --output ./results -profile docker
```

**What this default run does and costs.** It runs fastp, removes phiX and human reads with bwa-mem3, estimates sequencing coverage with Nonpareil, profiles the reads with sylph, assembles each sample with MEGAHIT and writes the results summary. On first use it downloads the human genome (940 MB) and the sylph database (19 GB), and builds the human bwa-mem3 index (~75 GB RAM, once per `--databases_dir`, by default `<output>/../databases`); decontamination and MEGAHIT request 72 GB each, so set `--max_memory` to the machine's RAM. The two usual changes: `--assembly_mode none` for read-level results only, and `--read_aligner bowtie2` on a machine with less than ~75 GB of RAM.

Profiles combine an engine with a site: `-profile docker`, `singularity`, `apptainer` or `podman`; `-profile slurm,singularity`; `-profile aws,docker` (all profiles: [manual, section 6.2](docs/manual.md#62-available-profiles); clusters and clouds: [docs/deployment.md](docs/deployment.md)). Every option: `nextflow run main.nf --help` and [docs/parameters.md](docs/parameters.md).

## Results

Start with `results/00_summary/`: `results_summary.html` opens the report (it works offline), `results_summary.pdf` is the same report typeset, `results_summary.xlsx` holds every table, and `README.md` there describes every table and figure. The numbered folders (`01_quality_control/` ... `07_functional_annotation/`) hold each stage's own outputs. What every folder and file holds, and how to read it: [manual, section 7](docs/manual.md#7-outputs-and-how-to-read-them). Several finished runs (batches) can be summarised together: [manual, section 9](docs/manual.md#9-merging-batches).

## Documentation

| Document | Contents |
|----------|----------|
| [`docs/manual.md`](docs/manual.md) | The user manual: requirements, installation and test, samplesheet, choosing what runs, running, outputs and how to read them, databases, merging batches, examples |
| [`docs/parameters.md`](docs/parameters.md) | Every parameter, with its type, default and options |
| [`docs/troubleshooting.md`](docs/troubleshooting.md) | Problems and their solutions |
| [`docs/deployment.md`](docs/deployment.md) | Running on SLURM and on AWS, Google Cloud and Azure |
| [`docs/WORKFLOW_DIAGRAM.md`](docs/WORKFLOW_DIAGRAM.md) | One diagram per stage, with the parameter that gates each step |
| [`docs/DECONTAMINATION_QUICK_REFERENCE.md`](docs/DECONTAMINATION_QUICK_REFERENCE.md) | Host and phiX removal: other and multiple host genomes, pre-built indexes |
| [`docs/RGI_WILDCARD_USAGE.md`](docs/RGI_WILDCARD_USAGE.md) | Preparing CARD and WildCARD by hand for RGI |
| [`docs/YAML_PARAMETERS_GUIDE.md`](docs/YAML_PARAMETERS_GUIDE.md) | Writing a `-params-file` YAML file: types, lists, and how it combines with the command line |
| [`docs/DISK_OPTIMIZATION.md`](docs/DISK_OPTIMIZATION.md) | Disk use: the `low_disk` profile, work-directory cleanup, publishing modes |
| [`docs/benchmarks/`](docs/benchmarks/) | Benchmarks behind the defaults (aligner, assembly, binning, multi-sample binning, contig taxonomy, report system) |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | For contributors: development setup, conventions, tests |

## Citation

If you use BugBuster, please cite the pipeline with its version (`nextflow run main.nf --version`) and repository, https://github.com/gene2dis/BugBuster, and cite the tools that ran: the Methods page of the run's report (`00_summary/report/`, and the PDF) lists every tool and database with its version and reference.

## Licence

BugBuster is released under the [MIT licence](LICENSE). The tools it runs and the databases it downloads keep their own licences and terms of use; for example, CARD (used by RGI and K-MARVEL) is free for academic use, and commercial use needs a licence from CARD.

## Credits

gene2dis/BUGBUSTER is developed by the Microbial Data Science Lab, Facultad de Ingeniería, Universidad del Desarrollo. Its development was led by Francisco A. Fuentes 

We thank the following people for their extensive assistance in the development of this pipeline:

- Francisco A. Fuentes
- Juan A. Ugalde
- Carolina Curiqueo

If you have any question of how to use the pipeline, you can contact the developer at the mail ffuentessantander@gmail.com. We will be happy to answer your questions!
