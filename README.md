![image](https://github.com/user-attachments/assets/a10c01f6-ef6c-40c4-a4ac-26a0c4f87564)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.10.0-23aa62.svg)](https://www.nextflow.io/)
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
2. Optionally remove samples that fall below a minimum read count after quality filtering (`--min_read_sample`, default `0` = filter disabled).
3. Remove host contamintant reads. [`Bowtie2`](https://github.com/BenLangmead/bowtie2)
4. If requested Antibiotic resistance prediction at read level using KARGA and KARGVA [`KARGA`](https://github.com/DataIntellSystLab/KARGA), [`KARGVA`](https://github.com/DataIntellSystLab/KARGVA)
5. If requested AMR gene prediction with pathogen-of-origin analysis using RGI [`RGI`](https://github.com/arpcard/rgi)
6. Normalization of predicted genes by estimating cell number with ARGs-OAP. [`ARGs-OAP`](https://github.com/xinehc/args_oap)
7. If requested read-level functional profiling with Woltka against the Web of Life (WoLr2) genomes: Bowtie2 alignment, ORF classification and KO / EC / COG / Pfam / MetaCyc pathway read counts with copies per genome equivalent, reported separately from the assembly-based tables (needs no assembly). [`Woltka`](https://github.com/qiyunzhu/woltka), [`Web of Life`](https://biocore.github.io/wol/)
   - Or, as the alternative read-level backend, SUPER-FOCUS against its DB_90 SEED subsystem database (DIAMOND or MMseqs2): SEED subsystem levels 1-3 read counts, reported separately as well (no copies per genome equivalent: a SEED hit has no gene length). [`SUPER-FOCUS`](https://github.com/metageni/SUPER-FOCUS)
   - Or HUMAnN 4 (alpha 4.0.0a2) with a MetaPhlAn 4.1.2 prescreen: MetaCyc pathway, KO and EC abundances in RPK with copies per genome equivalent, reported separately as well. [`HUMAnN`](https://github.com/biobakery/humann), [`MetaPhlAn`](https://github.com/biobakery/MetaPhlAn)
8. Taxonomic profile [`Kraken2`](https://ccb.jhu.edu/software/kraken2/) or [`Sourmash`](https://sourmash.readthedocs.io/en/latest/index.html)
9. Abundance estimation [`Bracken`](https://github.com/jenniferlu717/Bracken)
10. Read traceback and taxonomic reports (abundance tables and Phyloseq-compatible outputs).
11. Genome assembly [`Megahit`](https://github.com/voutcn/megahit)
12. Contig filter [`BBmap`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
13. Taxonomic annotation of contigs using Blastn and BlobTools. [`BlobTools`](https://github.com/DRL/blobtools), [`Blast`](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
14. ORF prediction in contigs with Pyrodigal (one shared gene-calling pass feeding DeepARG and functional annotation). [`Pyrodigal`](https://github.com/althonos/pyrodigal)
15. If requested contig-level functional annotation with eggNOG-mapper, CAZy annotation with run_dbcan, per-gene abundance quantification with featureCounts, average genome size estimation with MicrobeCensus, and study-level TPM and copies-per-genome-equivalent tables per functional ontology (KO, COG, EC, Pfam, CAZy — with the eggNOG and dbCAN CAZy calls reported separately). [`eggNOG-mapper`](https://github.com/eggnogdb/eggnog-mapper), [`run_dbcan`](https://github.com/bcb-unl/run_dbcan), [`featureCounts`](https://subread.sourceforge.net), [`MicrobeCensus`](https://github.com/snayfach/MicrobeCensus)
16. Prediction of resistance genes at the contig level with DeepARG. [`DeepARG`](https://github.com/gaarangoa/deeparg)
17. Contig reports, scatter plot of taxonomy at Phylum level and scatter plot of resistance genes in contigs.
18. Binning with user-selectable tools (default: SemiBin; options: [`Metabat2`](https://bitbucket.org/berkeleylab/metabat/src/master/), [`SemiBin`](https://github.com/BigDataBiology/SemiBin), [`COMEBin`](https://github.com/ziyewang/COMEBin))
19. Binning refinement with [`MetaWrap`](https://github.com/bxlab/metaWRAP) (only when ≥2 binners are selected)
20. Bin quality prediction [`CheckM2`](https://github.com/chklovski/CheckM2)
21. Bin taxonomic prediction [`GTDB-TK`](https://github.com/Ecogenomics/GTDBTk)
22. Bin reports.
23. If requested MAG-level functional annotation of the refined bins with Bakta, one annotation per bin (requires ≥2 binners so only MetaWRAP quality-filtered bins are annotated). [`Bakta`](https://github.com/oschwengers/bakta)
24. If requested ARG clustering [`mmseqs2`](https://github.com/soedinglab/MMseqs2)
25. Assembly modes: "coassembly", "assembly", "none"

## Quick Start

### Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=25.10.0`)
- Container runtime: [Docker](https://docs.docker.com/engine/installation/), [Singularity](https://sylabs.io/guides/), or [Podman](https://podman.io/)

### Basic Usage

```bash
# Clone the repository
git clone https://github.com/gene2dis/BugBuster.git
cd BugBuster

# Create your samplesheet from the bundled example and edit it
# with the paths to your reads (see "Samplesheet Format" below)
cp assets/samplesheet_example.csv samplesheet.csv

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
| `apptainer` | Run with Apptainer containers |
| `conda` | Run with Conda environments (containers are recommended; most modules only define containers) |
| `slurm` | Run on SLURM HPC cluster |
| `aws` | Run on AWS Batch |
| `gcp` | Run on Google Cloud |
| `azure` | Run on Azure Batch |
| `low_disk` | Deletes the work dir after a successful run, for limited disk space (completed runs are not resumable; interrupted ones are) |
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

See [docs/deployment.md](docs/deployment.md) for detailed deployment instructions, and [docs/troubleshooting.md](docs/troubleshooting.md) for solutions to common problems.

## Databases

All databases are automatically downloaded on first use and stored in `<output_dir>/../databases/` by default (configurable via `--databases_dir`), separate from the results directory. Only databases required for your selected pipeline options will be downloaded.

You can use custom databases by specifying paths with `--custom_*` parameters (see table below).

### Automatic Download Databases

| Database | Size | Used For | Trigger Parameter | Custom Path Parameter |
|----------|------|----------|-------------------|----------------------|
| **phiX174 Genome** | 5.4 kB | PhiX contamination removal (combined index built locally) | `quality_control=true` | `--custom_phiX_fasta` |
| **Human Host Genome (T2T-CHM13v2.0)** | 940 MB | Host read removal (combined index built locally) | `quality_control=true` | `--custom_host_fasta` (or `--custom_decontamination_index` for a pre-built index) |
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
| **GTDB-TK r232** | ~61 GB download | Bin taxonomic classification (GTDB-Tk 2.7.2; only R232 data works) | `include_binning=true` | `--custom_gtdbtk_db` |
| **eggNOG 7 (emapper-3.0)** | 44 GB | Contig functional annotation (eggNOG-mapper v3; requires the singularity/apptainer profile while v3 is in beta) | `contig_level_functional=true` | `--custom_eggnog_db` |
| **dbCAN (db_v5-2-9_5-5-2026)** | 7.4 GB | CAZy annotation of predicted proteins (run_dbcan v5) | `contig_level_functional=true` (and `functional_cazy=true`, the default) | `--custom_dbcan_db` |
| **Bakta DB v6.0** | full 31.9 GB / light 1.3 GB download | MAG (bin) annotation with Bakta; flavor via `--bakta_db` (light is explicit and recorded in provenance) | `mag_level_functional=true` | `--custom_bakta_db` |
| **NCBI COG 2024 definitions** | 410 KB | Maps eggNOG's and WoLr2's COG ids to COG functional categories | `contig_level_functional=true` or `read_level_functional=woltka` | `--custom_cog_db` |
| **Web of Life WoLr2** | ~94 GB | Read-level functional profiling with Woltka (alignment needs ≥ 68 GB RAM) | `read_level_functional='woltka'` | `--custom_woltka_db` |
| **SUPER-FOCUS DB_90** | ~0.74 GB (DIAMOND) / ~0.9 GB (MMseqs2) download, selected aligner only | Read-level SEED subsystem profiling with SUPER-FOCUS | `read_level_functional='superfocus'` | `--custom_superfocus_db` |
| **HUMAnN 4 databases + MetaPhlAn vOct22** | ~71 GB (full ChocoPhlAn 44.8 GB) or ~33 GB (EC-filtered ChocoPhlAn 6.9 GB), plus UniRef90 EC-filtered 0.94 GB, utility mapping 2.8 GB and the MetaPhlAn database with its Bowtie2 index ~22.5 GB | Read-level profiling with HUMAnN 4.0.0a2 | `read_level_functional='humann'` | `--custom_humann_db` |

### Database Sources

- **phiX_index**: [`phage phiX174 genome`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=genbank)
- **host_db**: [`T2T-CHM13v2.0 genome (chm13v2.0.fa.gz)`](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz)
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
- **eggnog_db**: [`emapper-3.0 data`](https://data.cgmlab.org/eggnog-mapper/emapper-3.0/data/)
- **dbcan_db**: [`dbCAN S3 release db_v5-2-9_5-5-2026`](https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026/)
- **woltka_db**: [`Web of Life release 2 (WoLr2)`](https://ftp.microbio.me/pub/wol2/)
- **superfocus_db**: [`SUPER-FOCUS DB_90 (figshare, CC0): DIAMOND`](https://doi.org/10.25451/flinders.25009748.v2) / [`MMseqs2`](https://doi.org/10.25451/flinders.25009751.v2)
- **humann_db**: [`HUMAnN v4_alpha databases`](https://github.com/biobakery/humann) (ChocoPhlAn, UniRef90 EC-filtered, utility mapping) + [`MetaPhlAn mpa_vOct22_CHOCOPhlAnSGB_202403`](https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/)
- **cog_db**: [`NCBI COG 2024 definitions (cog-24.def.tab)`](https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/)

## Samplesheet Format

Create a CSV file with your sample information (a template is provided at
`assets/samplesheet_example.csv`):

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
| `--contig_level_functional` | `false` | Contig functional annotation (Pyrodigal + eggNOG-mapper v3 + run_dbcan CAZy + featureCounts + TPM/CPGE tables); needs an assembly and `-profile singularity` or `apptainer` |
| `--mag_level_functional` | `false` | Bakta annotation of every refined bin (needs `--include_binning` and ≥2 `--binners`) |
| `--read_level_functional` | `none` | Read-level functional profiling: `woltka`, `superfocus`, `humann` or `none` (needs no assembly) |
| `--microbecensus` | `true` | MicrobeCensus average genome size for CPGE normalization (with a contig or read functional branch) |
| `--functional_cazy` | `true` | run_dbcan CAZy annotation in the contig functional branch |
| `--arg_bin_clustering` | `false` | Bin-level ARG prediction and clustering |

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
| `--bbmap_length` | `1000` | Minimum contig length after filtering |
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
├── pipeline_info/              # Execution reports, logs, contig filtering summary
├── clean_reads/                # Decontaminated reads (only if --store_clean_reads)
├── 01_quality_control/         # FastP reports (if quality_control=true) + read-count summary
│   ├── fastp/                  # Per-sample FastP reports
│   └── summary/                # Aggregated statistics
├── 02_taxonomy/                # Taxonomic profiling (if taxonomic_profiler != 'none')
│   ├── kraken2/ or sourmash/   # Per-sample profiler results
│   ├── bracken/                # Bracken abundance estimates (kraken2 only)
│   ├── tables/                 # Abundance tables
│   ├── phyloseq/               # Phyloseq objects (*.h5, *.RDS)
│   └── figures/                # Taxonomy plots
├── 03_assembly/                # Genome assembly (if assembly_mode != 'none')
│   ├── per_sample/             # Per-sample assemblies
│   └── coassembly/             # Co-assembly results
├── 04_binning/                 # Metagenomic binning (if include_binning=true)
│   ├── per_sample/{sample}/    # Assembly mode: raw_bins/, refined_bins/ (if >= 2 binners), quality/ (CheckM2), taxonomy/ (GTDB-TK)
│   ├── per_sample/coassembly/  # Co-assembly mode: quality/ (CheckM2) and taxonomy/ (GTDB-TK) only
│   ├── coassembly/             # Co-assembly mode: raw_bins/, refined_bins/ (if >= 2 binners), coverage/, summary/
│   ├── quality/summary/        # Aggregated quality reports (assembly mode)
│   └── taxonomy/summary/       # Aggregated taxonomy reports (assembly mode)
├── 05_arg_prediction/          # ARG predictions
│   ├── read_level/             # Read-level predictions
│   │   ├── karga/              # KARGA results (if read_arg_prediction=true)
│   │   ├── kargva/             # KARGVA results (if read_arg_prediction=true)
│   │   ├── args_oap/           # ARGs-OAP results (if read_arg_prediction=true)
│   │   ├── summary/            # Normalized ARG report (if read_arg_prediction=true)
│   │   ├── rgi/                # RGI per-sample results (if rgi_prediction=true)
│   │   ├── rgi_kmer/           # RGI pathogen-of-origin (if rgi_prediction=true)
│   │   └── rgi_summary/        # RGI aggregated reports (if rgi_prediction=true)
│   ├── contig_level/           # Contig-level ARG (if contig_tax_and_arg=true)
│   │   ├── deeparg/            # Per-sample DeepARG predictions (proteins from 07_functional_annotation/gene_calling/; coassembly/ under co-assembly)
│   │   ├── summary/            # Combined tax + ARG report
│   │   └── figures/            # ARG scatter plots
│   └── bin_level/              # Bin-level ARG (if arg_bin_clustering=true)
│       ├── proteins/           # Per-bin ORFs (Prodigal)
│       ├── deeparg/            # Per-bin DeepARG predictions
│       └── clustering/         # MMseqs2 ARG clusters
├── 06_contig_taxonomy/         # Contig taxonomy (if contig_tax_and_arg=true)
│   └── figures/                # BlobTools plots
└── 07_functional_annotation/   # Functional annotation
    ├── gene_calling/{sample}/  # Pyrodigal ORFs (if contig_tax_and_arg or contig_level_functional; coassembly/ under co-assembly)
    ├── eggnog/{sample}/        # eggNOG-mapper annotations (if contig_level_functional=true; coassembly/ under co-assembly)
    ├── dbcan/{sample}/         # run_dbcan CAZy calls + substrates (if contig_level_functional=true and functional_cazy=true; coassembly/ under co-assembly)
    ├── gene_abundance/{sample}/# featureCounts per-gene counts (if contig_level_functional=true; per sample in both assembly modes)
    ├── microbecensus/{sample}/ # MicrobeCensus average genome size (if microbecensus=true and contig_level_functional=true or read_level_functional != none; per sample in both assembly modes)
    ├── summary/                # Study-level TPM/CPGE tables + annotated fraction + AGS summary (if contig_level_functional=true); read_* tables (if read_level_functional != none)
    ├── mags/{sample}/          # Bakta per-bin MAG annotation (if mag_level_functional=true; needs binning with >= 2 binners; coassembly/ under co-assembly)
    ├── reads/woltka/{sample}/  # Woltka read-level ORF and function tables (if read_level_functional=woltka)
    ├── reads/superfocus/{sample}/ # SUPER-FOCUS SEED level 1-3 tables + raw outputs (if read_level_functional=superfocus)
    └── reads/humann/{sample}/  # HUMAnN gene families / reactions / pathways (RPK), MetaPhlAn profile, KO/EC/MetaCyc tables (if read_level_functional=humann)
```

**Database Storage**: Databases are stored separately at `<output_dir>/../databases/` by default (configurable via `--databases_dir`).

For detailed output descriptions, see [`docs/manual.md`](docs/manual.md#8-output-structure).

## Documentation

| Document | Contents |
|----------|----------|
| [`docs/manual.md`](docs/manual.md) | Complete user manual: installation, databases, parameters, examples, output structure |
| [`docs/parameters.md`](docs/parameters.md) | Full parameter reference with types, defaults, and examples |
| [`docs/troubleshooting.md`](docs/troubleshooting.md) | Solutions to common problems |
| [`docs/deployment.md`](docs/deployment.md) | Deployment on HPC (SLURM) and cloud (AWS, GCP, Azure) |
| [`docs/WORKFLOW_DIAGRAM.md`](docs/WORKFLOW_DIAGRAM.md) | Pipeline stage diagram and the parameters that gate each stage |
| [`docs/DECONTAMINATION_QUICK_REFERENCE.md`](docs/DECONTAMINATION_QUICK_REFERENCE.md) | Host/phiX decontamination options and custom reference genomes |
| [`docs/RGI_WILDCARD_USAGE.md`](docs/RGI_WILDCARD_USAGE.md) | Manual CARD/WildCARD database preparation for RGI |
| [`docs/YAML_PARAMETERS_GUIDE.md`](docs/YAML_PARAMETERS_GUIDE.md) | Running the pipeline with `-params-file` YAML files |
| [`docs/DISK_OPTIMIZATION.md`](docs/DISK_OPTIMIZATION.md) | Reducing disk usage (`low_disk` profile, work-dir cleanup) |

## Credits

gene2dis/BUGBUSTER was originally written by the Microbial Data Science Lab, Center for Bioinformatics and Integrative Biology, Universidad Andres Bello. Its development was led by Francisco A. Fuentes 

We thank the following people for their extensive assistance in the development of this pipeline:

- Francisco A. Fuentes
- Juan A. Ugalde
- Carolina Curiqueo

If you have any question of how to use the pipeline, you can contact the developer at the mail ffuentessantander@gmail.com. We will be happy to answer your questions!
