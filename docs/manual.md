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
- **Functional Annotation**: Contig-level with eggNOG-mapper v3 (Pyrodigal gene calling + eggNOG 7 orthology transfer), run_dbcan v5 CAZy annotation with substrate prediction, per-gene abundance quantification (featureCounts), average genome size estimation (MicrobeCensus) and study-level TPM and copies-per-genome-equivalent tables per functional ontology (KO, COG, EC, Pfam, CAZy; requires singularity/apptainer, see the note under the feature toggles), or with MetaCerberus; MAG-level with Bakta (per-bin annotation of MetaWRAP-refined bins); optional read-level profiling with Woltka against the Web of Life (WoLr2) genomes (KO, EC, COG, Pfam, MetaCyc pathway read counts and copies per genome equivalent) or with SUPER-FOCUS against its SEED subsystem database (SEED subsystem levels 1–3 read counts), one backend per run, reported separately from the contig branch; needs no assembly

### Pipeline Workflow

```
Reads → QC (FastP) → Host Removal (Bowtie2) → Taxonomy (Kraken2/Sourmash)
                                            ├──→ Read-level functions (Woltka/WoLr2 or SUPER-FOCUS/SEED, optional)
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
                              ↓
                      MAG Annotation (Bakta, optional)
```

---

## 2. Requirements

### Software Requirements

| Software | Minimum Version | Purpose |
|----------|-----------------|---------|
| Nextflow | ≥25.10.0 | Workflow engine |
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

BugBuster automatically downloads required databases on first use. Databases are stored in `<output_dir>/../databases/` by default (configurable via `--databases_dir`), separate from the results directory, so they can be reused across pipeline runs.

### Automatic Download Databases

| Database | Size | Used For | Download Trigger |
|----------|------|----------|------------------|
| **phiX174 genome (FASTA)** | 5.4 kB | PhiX contamination removal | `quality_control=true` |
| **Human host genome T2T-CHM13v2.0 (FASTA)** | 940 MB | Host read removal | `quality_control=true` |
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
| **GTDB-TK r232** | ~61 GB download | Bin taxonomic classification (GTDB-Tk 2.7.2; only R232 data works) | `include_binning=true` |
| **eggNOG 7 (emapper-3.0)** | 44 GB | Contig functional annotation | `contig_level_functional=true` |
| **dbCAN (db_v5-2-9_5-5-2026)** | 7.4 GB | CAZy annotation of predicted proteins | `contig_level_functional=true` (unless `functional_cazy=false`) |
| **NCBI COG 2024 definitions (cog-24.def.tab)** | 410 KB | Maps eggNOG's COG ids to COG functional categories | `contig_level_functional=true` |
| **Bakta DB v6.0 full** | 31.9 GB download | MAG (bin) annotation | `mag_level_functional=true` |
| **Bakta DB v6.0 light** | 1.3 GB download | MAG (bin) annotation, reduced annotation sources | `mag_level_functional=true` with `bakta_db='v6.0-light'` (explicit choice, recorded in provenance) |
| **Web of Life WoLr2** | ~94 GB (Bowtie2 index 93.6 GB + coordinates/maps ~0.6 GB) | Read-level functional profiling (Woltka); alignment needs ≥ 68 GB RAM | `read_level_functional='woltka'` |
| **SUPER-FOCUS DB_90** | ~0.74 GB download (~1.9 GB unpacked) for DIAMOND, or ~0.9 GB (~2.5 GB unpacked) for MMseqs2; only the selected aligner's archive | Read-level SEED subsystem profiling (SUPER-FOCUS) | `read_level_functional='superfocus'` |

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
wget -O /shared/databases/bugbuster/gtdbtk_r232.tar.gz \
    https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
tar -xzf /shared/databases/bugbuster/gtdbtk_r232.tar.gz -C /shared/databases/bugbuster/

# Download the Bakta database v6.0 (full, 31.9 GB; use db-light.tar.xz for
# the 1.3 GB light DB — note the light DB changes annotation results)
wget -P /shared/databases/bugbuster/bakta/ \
    https://zenodo.org/record/14916843/files/db.tar.xz
tar -xJf /shared/databases/bugbuster/bakta/db.tar.xz -C /shared/databases/bugbuster/bakta/

# Download the Web of Life (WoLr2) subset the Woltka read-level backend uses
# (~94 GB), keeping the FTP layout so the directory works as --custom_woltka_db.
# The .md5 files are checksums of the UNCOMPRESSED content.
BASE=https://ftp.microbio.me/pub/wol2
DEST=/shared/databases/bugbuster/wol2
get() { mkdir -p "$DEST/$(dirname "$1")"; wget -c --tries=10 -nv -O "$DEST/$1" "$BASE/$1"; }
for f in WoLr2.1.bt2l WoLr2.2.bt2l WoLr2.3.bt2l WoLr2.4.bt2l WoLr2.rev.1.bt2l WoLr2.rev.2.bt2l; do
    get databases/bowtie2/$f; done
for f in coords.txt length.map; do get proteins/$f.xz; get proteins/$f.md5; done
for f in orf-to-ko.map.xz orf-to-ko.map.md5 ko-to-ec.map ko-to-cog.map ko_name.txt; do get function/kegg/$f; done
for f in orf-to-protein.map.xz orf-to-protein.map.md5 protein-to-enzrxn.map enzrxn-to-reaction.map \
         reaction-to-pathway.map pathway_name.txt; do get function/metacyc/$f; done
for f in orf-to-pfam.map.xz orf-to-pfam.map.md5 pfam_name.txt; do get function/pfam/$f; done
for x in proteins/coords.txt proteins/length.map function/kegg/orf-to-ko.map \
         function/metacyc/orf-to-protein.map function/pfam/orf-to-pfam.map; do
    [ "$(xz -dc $DEST/$x.xz | md5sum | cut -d' ' -f1)" = "$(cut -d' ' -f1 $DEST/$x.md5)" ] \
        && echo "OK  $x" || echo "BAD $x"; done

# SUPER-FOCUS DB_90 for the SUPER-FOCUS read-level backend (figshare, CC0).
# The database ROOT is the directory CONTAINING db/; pass it as
# --custom_superfocus_db. DIAMOND archive shown (~0.74 GB); for
# --superfocus_aligner mmseqs2 use file 44075237 into db/static/mmseqs2 instead.
mkdir -p /shared/databases/bugbuster/superfocus/db/static/diamond
wget -O db90.zip https://ndownloader.figshare.com/files/44075225
unzip -j db90.zip -d /shared/databases/bugbuster/superfocus/db/static/diamond
wget -O /shared/databases/bugbuster/superfocus/db/database_PKs.txt \
    https://raw.githubusercontent.com/metageni/SUPER-FOCUS/739404db8816de967cd4ac3e0d9effcdab7f1489/superfocus_app/db/database_PKs.txt

# NCBI COG 2024 definitions table (~410 KB; maps eggNOG's COG ids to COG
# functional categories in the contig branch)
wget -P /shared/databases/bugbuster/cog/ \
    https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.def.tab

# Download human host genome (T2T-CHM13v2.0); the pipeline builds the
# combined phiX + host Bowtie2 index from FASTA on first use
wget -P /shared/databases/bugbuster/ \
    https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
```

### Using Custom Databases

Specify custom database paths to skip automatic downloads:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --custom_kraken_db /shared/databases/bugbuster/kraken2_standard8 \
    --custom_host_fasta /shared/databases/bugbuster/chm13v2.0.fa.gz \
    --custom_checkm2_db /shared/databases/bugbuster/checkm2/uniref100.KO.1.dmnd \
    --custom_gtdbtk_db /shared/databases/bugbuster/release232 \
    --custom_bakta_db /shared/databases/bugbuster/bakta/db \
    --custom_woltka_db /shared/databases/bugbuster/wol2 \
    --custom_superfocus_db /shared/databases/bugbuster/superfocus \
    --custom_cog_db /shared/databases/bugbuster/cog/cog-24.def.tab \
    -profile docker
```

### Using Custom Host Genomes

To filter reads from non-human hosts, point the pipeline at the host genome
FASTA (plain or gzipped) — it builds the combined phiX + host Bowtie2 index
automatically:

```bash
# Download your host genome
wget -O host_genome.fasta.gz <URL_TO_HOST_GENOME>

# Use in pipeline
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --custom_host_fasta ./host_genome.fasta.gz \
    -profile docker
```

If you already have a pre-built combined (phiX + host) Bowtie2 index, pass its
directory with `--custom_decontamination_index` instead.

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
| `sample` | Yes | Unique sample identifier: letters, digits, underscore, dot or hyphen, starting with a letter or digit (it becomes file/directory names) |
| `r1` | Yes | Absolute path to forward reads (R1) in FASTQ/FASTQ.GZ format |
| `r2` | Yes | Absolute path to reverse reads (R2) in FASTQ/FASTQ.GZ format |
| `s` | No | Absolute path to singleton reads (optional, leave empty if none) |

Singleton reads are processed in both QC modes: with `--quality_control true`
(default) they are trimmed by a dedicated single-end fastp run, decontaminated
alongside the paired reads, and carried into downstream assembly/profiling
steps that accept them. A few caveats: the `--min_read_sample` threshold counts
paired reads only; Kraken2 profiles the paired reads only (it runs in
`--paired` mode); and read-level RGI ARG prediction (`--rgi_prediction`) also
uses the paired reads only (`rgi bwt` accepts a single read pair). Sourmash
and assembly use singletons as well.

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

Parameters are validated against `nextflow_schema.json` at startup (nf-schema). Passing a parameter that the pipeline does not declare — including misspelled or removed ones — aborts the run with an "unrecognised parameter" error. See [`parameters.md`](parameters.md) for the complete reference.

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
| `--contig_level_functional` | `false` | `true`, `false` | Contig functional annotation (Pyrodigal + eggNOG-mapper v3 + run_dbcan CAZy + featureCounts gene abundance + TPM/CPGE summary tables); see engine note below |
| `--microbecensus` | `true` | `true`, `false` | MicrobeCensus average genome size for CPGE normalization (with the contig or read functional branch; failure is non-fatal — affected samples get no CPGE: TPM-only in the contig tables, native read counts only in the read tables) |
| `--functional_cazy` | `true` | `true`, `false` | run_dbcan v5 CAZy annotation of the predicted proteins (only with the functional branch; reported separately from the eggNOG CAZy calls, never merged) |
| `--dbcan_consensus` | `recommended` | `recommended`, `any` | Which dbCAN calls feed the summary tables (`recommended` = supported by ≥2 tools; `any` = per-tool union) |
| `--mag_level_functional` | `false` | `true`, `false` | MAG-level functional annotation: Bakta on every refined bin (requires `--include_binning` and ≥2 `--binners`; see the note below) |
| `--read_level_functional` | `none` | `woltka`, `superfocus`, `none` | Read-level functional profiling backend (needs no assembly; see the notes below) |
| `--woltka_uniq` | `false` | `true`, `false` | Woltka: leave reads whose reported alignments hit several ORFs unassigned instead of dividing them 1/k |
| `--superfocus_aligner` | `diamond` | `diamond`, `mmseqs2` | SUPER-FOCUS search backend (DIAMOND blastx or MMseqs2); selects which DB_90 archive is downloaded |
| `--contig_level_metacerberus` | `false` | `true`, `false` | Functional annotation with MetaCerberus |
| `--arg_bin_clustering` | `false` | `true`, `false` | ARG clustering for HGT inference |
| `--min_read_sample` | `0` | Integer ≥ 0 | Minimum reads required after QC |

> **Engine note for `--contig_level_functional`:** eggNOG-mapper v3 is in beta and its
> authors publish only an Apptainer image, so this branch requires
> `-profile singularity` or `-profile apptainer`. Nextflow uses one container engine
> per run, which means enabling this branch switches the **entire run** to that
> engine — every other tool runs from the same pinned images, automatically
> converted, so results are unchanged (expect a one-time image-conversion delay on
> first run). Because the cloud profiles (`aws`, `gcp`, `azure`) are docker-based,
> this branch cannot run on cloud batch executors during the beta. Runs without
> this flag are unaffected on every engine. Docker/cloud support returns when
> eggNOG-mapper v3.0.0 final is released on bioconda.

> **Note for `--mag_level_functional`:** Bakta annotates the refined bins, one task
> per bin, publishing per-bin GFF3/GBFF/FAA/FNA/TSV/summary files under
> `07_functional_annotation/mags/<sample>/`. It requires `--include_binning` **and at
> least two `--binners`**: MetaWRAP refinement and its completeness/contamination
> quality filter (defaults 50/10) only run when ≥2 binners are selected, and only
> quality-filtered bins are annotated. Bins are expected to be **bacterial** — Bakta
> is a bacterial annotator, and the pipeline does not detect or exclude
> archaeal/eukaryotic/viral bins; interpret annotations of such bins with caution
> (check the GTDB-Tk bin taxonomy report). This branch is independent of
> `--contig_level_functional` and, unlike it, runs on any container engine
> (docker included).

> **Note for `--read_level_functional woltka`:** host-removed reads are aligned with
> Bowtie2 to the Web of Life release 2 (WoLr2, 15,953 genomes) using the SHOGUN
> multi-hit settings (`-k 16`, the WoL/Qiita standard), and Woltka assigns each read
> to the WoLr2 ORF it overlaps (≥ 80 % of the read inside the ORF). Mates are counted
> as separate reads, and a read whose reported alignments (up to 16 within the
> score threshold) hit k ORFs contributes 1/k to each
> (`--woltka_uniq` leaves such reads unassigned instead). ORF counts are then
> summarized to KO, EC (via KO), COG ortholog groups (via KO; e.g. `COG0604` — not
> the single-letter COG categories of the contig branch), Pfam and MetaCyc pathways
> from the WoLr2 maps. **A read contributes its full weight once to each distinct
> term its ORF carries** (an ORF with two KOs counts toward both; the same intentional
> double counting as the contig branch). The branch runs on any container engine,
> needs no assembly (it also works with `--assembly_mode none`), and its tables are
> reported **separately** (`read_*` files). Woltka identifies reads by name, so input
> reads must have unique names (normal for sequencer output; some simulated test
> datasets reuse names, and Woltka then counts same-named reads once). Read-level
> profiling recovers the unassembled fraction but over-predicts, so it is never
> merged with the assembly-based tables. The WoLr2 database is ~94 GB and
> alignment needs **≥ 68 GB RAM** per task; pre-download it once and pass
> `--custom_woltka_db` (Section 4, Manual Database Download).

> **Note for `--read_level_functional superfocus`:** each sample's host-removed R1, R2
> and singleton reads are concatenated into one query (mates counted as separate
> reads) and profiled with SUPER-FOCUS 1.8 against its DB_90 SEED subsystem cluster
> database, using DIAMOND blastx (default) or MMseqs2 (`--superfocus_aligner
> mmseqs2`). For each read, its equal-best-e-value hits passing the SUPER-FOCUS
> defaults (≥ 60 % identity, ≥ 15 aa alignment, e-value ≤ 1e-5; set in the
> `SUPERFOCUS` `ext.args`) are counted, and **each read with a hit contributes
> exactly 1**, divided 1/k across the k distinct SEED (subsystem, function)
> assignments of its best hits. The pipeline sums these function-level counts to SEED
> subsystem levels 1–3 (ontologies `seed_level1`, `seed_level2`, `seed_level3`;
> accessions are path-qualified — `L1`, `L1 | L2`, `L1 | L2 | L3` — because the
> level-2 placeholder `-` occurs under many level-1 categories). SEED is never
> mapped to KO or EC, and the read tables are never merged with the contig branch.
> **No copies per genome equivalent are produced:** a SEED hit carries no gene
> length, so no RPK (and no CPGE) exists — `abundance_cpge` stays blank and
> `cpge_status` is `not_applicable`; MicrobeCensus still runs and its AGS / genome
> equivalents are reported in `read_sample_summary.tsv`. The cluster level is fixed
> at DB_90 (~0.74 GB DIAMOND / ~0.9 GB MMseqs2 download, only the selected
> aligner's archive). MMseqs2 builds its k-mer index for the database at every run
> (~13 GB RAM, ~45 s on 16 CPUs), and in its fast mode its result for borderline
> reads can vary slightly with the thread count; DIAMOND results are stable. In its
> fast mode MMseqs2 was also somewhat less sensitive in a single spot check (one
> 2,000-read test sample: DIAMOND 1,011 reads hit, MMseqs2 924). The branch runs on
> any container engine (docker included) and needs no assembly.

### 6.3 Database Selection Options

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `--phiX_index` | `phiX174` | `phiX174` | PhiX genome for decontamination |
| `--host_db` | `human` | `human` | Host genome for read removal |
| `--kraken2_db` | `standard-8` | `standard-8`, `gtdb_220` | Kraken2 database version |
| `--sourmash_db` | `gtdb_220_k31` | `gtdb_220_k31` | Sourmash database version |
| `--checkm2_db` | `v3` | `v3` | CheckM2 database version |
| `--gtdbtk_db` | `release_232` | `release_232` | GTDB-TK database release (only R232 works with the pinned GTDB-Tk 2.7.2) |
| `--eggnog_db` | `emapper-3.0` | `emapper-3.0` | eggNOG 7 data for eggNOG-mapper v3 |
| `--dbcan_db` | `db_v5-2-9_5-5-2026` | `db_v5-2-9_5-5-2026` | dbCAN database release for run_dbcan v5 |
| `--bakta_db` | `v6.0-full` | `v6.0-full`, `v6.0-light` | Bakta database flavor; the light DB changes annotation results, so selecting it is always explicit and is recorded in provenance (`software_versions.yml`) |
| `--woltka_db` | `wolr2` | `wolr2` | Web of Life release for the Woltka read-level backend |
| `--superfocus_db` | `db90` | `db90` | SUPER-FOCUS database for the SUPER-FOCUS read-level backend (DB_90, figshare CC0) |
| `--cog_db` | `cog-24` | `cog-24` | NCBI COG definitions table mapping eggNOG's COG ids to COG functional categories |

### 6.4 Custom Database Paths

| Parameter | Description |
|-----------|-------------|
| `--custom_decontamination_index` | Path to a pre-built combined Bowtie2 index directory (host + phiX) |
| `--custom_phiX_fasta` | Path to a custom phiX genome FASTA (index is built by the pipeline) |
| `--custom_host_fasta` | Path to a custom host genome FASTA (index is built by the pipeline) |
| `--custom_kraken_db` | Path to custom Kraken2 database directory |
| `--custom_sourmash_db` | List of paths: `["kmer.zip", "lineages.csv"]` |
| `--custom_checkm2_db` | Path to CheckM2 database file |
| `--custom_gtdbtk_db` | Path to the directory directly containing the unarchived GTDB-Tk reference data (e.g. the extracted `release232/`); only R232 data works with the pinned GTDB-Tk 2.7.2, R220/R226 does not (see `docs/parameters.md`) |
| `--custom_deeparg_db` | Path to DeepARG database directory |
| `--custom_blast_db` | Path to BLAST NT database directory |
| `--custom_taxdump_files` | Path to NCBI taxdump directory |
| `--custom_karga_db` | Path to KARGA database FASTA file |
| `--custom_kargva_db` | Path to KARGVA database FASTA file |
| `--custom_rgi_card_db` | Path to pre-prepared CARD database directory |
| `--custom_rgi_wildcard` | Path to WildCARD directory (use with `--custom_rgi_card_db`) |
| `--custom_eggnog_db` | Path to eggNOG 7 data directory (emapper-3.0 layout, see `docs/parameters.md`) |
| `--custom_dbcan_db` | Path to dbCAN database directory (run_dbcan v5 layout, see `docs/parameters.md`) |
| `--custom_bakta_db` | Path to Bakta database directory (schema 6 layout, see `docs/parameters.md`) |
| `--custom_woltka_db` | Path to a local WoLr2 mirror in the FTP layout (see `docs/parameters.md`) |
| `--custom_superfocus_db` | Path to a SUPER-FOCUS database ROOT, the directory containing `db/` (see `docs/parameters.md`) |
| `--custom_cog_db` | Path to a local NCBI COG definitions table (`cog-24.def.tab` layout) |

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
| `--bracken_read_len` | `150` | Read length for Bracken estimation |
| `--bracken_tax_level` | `S` | Taxonomic level: D, P, C, O, F, G, S |
| `--sourmash_tax_rank` | `species` | Rank: `genus`, `species`, `strain` |

### 6.7 Assembly Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bbmap_length` | `1000` | Minimum contig length after BBMap filtering |

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

**Note**: For detailed manual CARD database preparation instructions, see [`docs/RGI_WILDCARD_USAGE.md`](RGI_WILDCARD_USAGE.md) ("Manual CARD Database Preparation").

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

### 6.12 MetaCerberus Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metacerberus_hmm` | `"KOFam_all, COG, VOG, PHROG, CAZy"` | HMM databases to use |
| `--metacerberus_minscore` | `25` | Minimum HMM score |
| `--metacerberus_evalue` | `1e-09` | Maximum E-value |

**Available HMM Databases:** `KOFam_all`, `KOFam_eukaryote`, `KOFam_prokaryote`, `COG`, `VOG`, `PHROG`, `CAZy`

**Note:** the value is passed verbatim to MetaCerberus's `--hmm` flag, so it must keep embedded double quotes to stay a single argument. On the command line, wrap it in single quotes: `--metacerberus_hmm '"KOFam_prokaryote, COG, CAZy"'`.

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

By default, databases are stored in a `databases/` directory at the same level as your output directory. This allows database reuse across multiple pipeline runs.

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
| `low_disk` | Minimize disk usage: automatic work-dir cleanup (`cleanup = true`, disables `-resume`) plus `--store_clean_reads` — see [`DISK_OPTIMIZATION.md`](DISK_OPTIMIZATION.md) |
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
│   ├── pipeline_dag_*.svg                      # Pipeline DAG
│   ├── software_versions.yml                   # Versions of every tool used in the run
│   └── contig_filtering_summary.txt            # Contig filtering summary (if assembly_mode != 'none')
├── clean_reads/                                # Decontaminated reads (only if --store_clean_reads)
│   └── {sample}/                               # Per-sample clean R1/R2/Singleton FASTQs
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
│   ├── bracken/                                # Bracken abundance estimates (kraken2 only)
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
│   │       ├── {sample}_filtered_contigs.fa    # Filtered contigs (≥ bbmap_length bp)
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
│   │   ├── deeparg/                            # DeepARG predictions per sample
│   │   │   └── {sample}/
│   │   │       └── *_contigs_deep_arg.out.mapping.ARG
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
└── 07_functional_annotation/                   # Functional annotation
    ├── gene_calling/                           # Pyrodigal ORF predictions on contigs
    │   └── {sample}/                           # (if contig_tax_and_arg=true or contig_level_functional=true)
    │       ├── {sample}.faa.gz
    │       ├── {sample}.fna.gz
    │       ├── {sample}.gff.gz
    │       └── {sample}.score.gz
    ├── eggnog/                                 # eggNOG-mapper v3 functional annotation
    │   └── {sample}/                           # (if contig_level_functional=true; needs singularity/apptainer)
    │       ├── {sample}.emapper.seed_orthologs
    │       └── {sample}.emapper.annotations
    ├── dbcan/                                  # run_dbcan v5 CAZy annotation
    │   └── {sample}/                           # (if contig_level_functional=true and functional_cazy=true)
    │       ├── {sample}.overview.tsv           # per-gene CAZy calls with per-tool columns retained
    │       ├── {sample}.dbCANsub_hmm_results.tsv  # dbCAN-sub results incl. substrate predictions
    │       ├── {sample}.dbCAN_hmm_results.tsv  # raw dbCAN HMM results (provenance)
    │       └── {sample}.diamond.out            # raw DIAMOND-vs-CAZy results (provenance)
    ├── gene_abundance/                         # featureCounts per-gene read counts
    │   └── {sample}/                           # (if contig_level_functional=true; per sample in both assembly modes)
    │       ├── {sample}.featureCounts.txt      # Geneid, coordinates, Length, read count
    │       └── {sample}.featureCounts.txt.summary  # assigned vs unassigned alignments
    ├── microbecensus/                          # MicrobeCensus average genome size
    │   └── {sample}/                           # (if microbecensus=true and contig_level_functional=true or read_level_functional != none)
    │       ├── {sample}.ags.tsv                # AGS, genome equivalents, total bases
    │       └── {sample}.microbecensus.txt      # raw MicrobeCensus output (provenance)
    ├── summary/                                # study-level tables (if contig_level_functional=true)
    │   ├── gene_annotations.tsv                # long format: one row per gene per functional term
    │   ├── gene_abundance.tsv                  # per-gene counts + TPM + CPGE per sample
    │   ├── function_abundance.tsv              # per-ontology TPM and CPGE (ko, cog, ec, pfam, cazy; backend column separates eggnog-mapper and run_dbcan rows)
    │   ├── function_wide_{ontology}_tpm.tsv    # wide TPM matrix per ontology (rows terms, columns samples; ko/cog/ec/pfam/cazy from eggNOG, cazy_dbcan from run_dbcan)
    │   ├── function_wide_{ontology}_cpge.tsv   # wide CPGE matrix per ontology (blank columns for samples without AGS)
    │   ├── annotated_fraction.tsv              # per-sample annotated fraction by count and abundance
    │   ├── ags_and_ge.tsv                      # per-sample AGS summary with status ok / unavailable
    │   ├── read_function_abundance.tsv         # read branch (if read_level_functional != none): same schema, source=reads
    │   ├── read_function_wide_{ontology}_native.tsv  # wide read-count matrix (woltka: ko, ec, cog, pfam, metacyc; superfocus: seed_level1/2/3)
    │   ├── read_function_wide_{ontology}_cpge.tsv    # wide CPGE matrix (woltka only; blank columns for samples without AGS)
    │   ├── read_annotated_fraction.tsv         # share of reads carrying a term, per ontology (woltka: of ORF-assigned reads; superfocus: of all reads)
    │   └── read_sample_summary.tsv             # reads assigned / ambiguous, AGS, CPGE status, tool + DB version
    ├── reads/                                  # read-level functional profiling, one backend per run
    │   ├── woltka/{sample}/                    # (if read_level_functional=woltka; needs no assembly)
    │   │   ├── {sample}.bowtie2.log            # alignment summary against WoLr2
    │   │   ├── {sample}.woltka_orf.tsv         # Woltka per-ORF read counts
    │   │   ├── {sample}.woltka_functions.tsv   # per-function read counts and RPK
    │   │   ├── {sample}.woltka_summary.tsv     # reads assigned / annotated per ontology
    │   │   └── {sample}.woltka_unassigned.tsv  # reads left unassigned by --woltka_uniq
    │   └── superfocus/{sample}/                # (if read_level_functional=superfocus; needs no assembly)
    │       ├── {sample}.superfocus_functions.tsv   # SEED level 1-3 read counts (seed_level1/2/3)
    │       ├── {sample}.superfocus_summary.tsv     # reads given to SUPER-FOCUS / with an accepted hit
    │       ├── {sample}.superfocus_all_levels_and_function.xls  # raw SUPER-FOCUS table (tab-separated)
    │       ├── {sample}.superfocus_subsystem_level_{1,2,3}.xls  # raw SUPER-FOCUS per-level tables (not used, see note)
    │       └── {sample}.superfocus.log         # SUPER-FOCUS run log
    ├── mags/                                   # Bakta MAG-level annotation, one file set per bin
    │   └── {sample}/                           # (if mag_level_functional=true; requires binning with >= 2 binners)
    │       ├── {sample}_{bin}.gff3             # annotation in GFF3
    │       ├── {sample}_{bin}.gbff             # annotation in GenBank flat file
    │       ├── {sample}_{bin}.faa              # protein sequences
    │       ├── {sample}_{bin}.fna              # replicon/contig sequences
    │       ├── {sample}_{bin}.tsv              # per-feature annotation table
    │       ├── {sample}_{bin}.txt              # per-bin annotation summary
    │       └── {sample}_{bin}.hypotheticals.tsv  # hypothetical-protein table (+ .hypotheticals.faa)
    └── contigs/                                # MetaCerberus contig-level annotation
        └── {sample}/                           # (if contig_level_metacerberus=true)
            └── {sample}_annotation_results/
```

> **Gene identifiers and assembly mode:** in the gene abundance tables, `Geneid`
> is `<contig>_<n>`, matching the protein ids in `gene_calling/` and the query
> ids in the eggNOG annotations. Under `--assembly_mode coassembly` every
> sample is counted against the same co-assembly gene set, so gene ids are
> shared across samples and gene-level comparisons between samples within the
> run are valid. Under `--assembly_mode assembly` each sample has its own
> assembly and gene set: gene ids are sample-specific, and only function-level
> results (not per-gene rows) may be compared across samples. Counts are
> read-level (each mate counted separately), not fragment-level; the
> multi-mapping policy is set by `--featurecounts_multimap`.

> **TPM normalization (`summary/` tables):**
>
> ```
> RPK_i = count_i / (length_i / 1000)
> TPM_i = RPK_i / (sum over all j of RPK_j) * 1e6
> ```
>
> The share of the sample's functional pool attributable to that gene.
> Compositional. Use for within-sample composition and for compositionally
> aware differential testing. TPM sums to 1e6 per sample; a sample whose genes
> attracted no reads at all has TPM 0 for every gene instead. In
> `function_abundance.tsv`, a gene carrying two terms of the same ontology
> contributes its full abundance to each — this double-counting is
> intentional, so ontology-level TPM totals can exceed 1e6. The `description`
> column is empty for eggNOG terms (eggNOG-mapper v3 dropped the Description
> field). Pfam terms are Pfam family names (e.g. `BPD_transp_1`): eggNOG-mapper
> v3 reports one value per domain hit with its coordinates appended
> (`BPD_transp_1_210_403`), which the aggregation strips, so a gene with a
> repeated domain counts once toward that family. COG terms are COG
> functional categories (one letter, e.g. `P` = inorganic ion transport and
> metabolism): eggNOG 7 reports most genes with a COG ortholog id instead of a
> category (e.g. `COG1629`), which is mapped to its category letters with
> NCBI's COG definitions table (`cog-24.def.tab`, downloaded with the branch);
> a COG with several categories contributes to each. The table's name is
> appended to `db_version` on the COG rows of `gene_annotations.tsv`.

> **Two CAZy backends (`--functional_cazy`, on by default):** CAZy calls come
> from both eggNOG-mapper (coarse, orthology-transferred) and run_dbcan v5
> (dedicated CAZyme annotation with family/subfamily resolution). They are
> reported **separately, never merged**: in the long tables the `db`
> (`eggnog_cazy` vs `dbcan`) and `backend` (`eggnog-mapper` vs `run_dbcan`)
> columns distinguish them, and the wide matrices are split into
> `function_wide_cazy_*` (eggNOG) and `function_wide_cazy_dbcan_*`
> (run_dbcan) — a merged matrix would double-count genes called by both
> tools. For CAZy-focused analyses prefer the dbCAN matrices. Which dbCAN
> calls feed the tables is set by `--dbcan_consensus` (`recommended` =
> supported by ≥2 of DIAMOND / dbCAN HMM / dbCAN-sub, the default; `any` =
> per-tool union); the published `overview.tsv` always retains the per-tool
> columns. In `annotated_fraction.tsv` the `cazy` rows count a CAZy call
> from either backend and the `cazy_dbcan` rows count run_dbcan alone. With
> `--functional_cazy false` the `cazy_dbcan` outputs are still written but
> empty. dbCAN rows carry no per-call e-value/score (the overview has no
> single per-call value), and their `description` column is empty. dbCAN
> accessions are kept exactly as the tool emits them: plain families
> (`GH13`), CAZy subfamilies (`GH5_4`), and dbCAN-sub subfamily cluster ids
> (`GH78_e118`) all occur.

> **Copies per genome equivalent (`cpge` / `abundance_cpge` columns and the
> `function_wide_*_cpge.tsv` matrices):**
>
> ```
> genome_equivalents = total_bases_sampled / average_genome_size_bp   (from MicrobeCensus)
> CPGE_i = RPK_i / genome_equivalents
> ```
>
> Approximate average copies of that gene per community member. Not
> compositional. Use when the question is whether an average cell carries the
> function, and for comparing across communities of different composition.
> TPM and CPGE answer different questions: they are not interchangeable and
> neither replaces the other — the most common error in this kind of analysis
> is treating one as the other.
>
> The genome-equivalents estimate comes from MicrobeCensus run on the
> host-removed reads (on by default with the branch; disable with
> `--microbecensus false`). MicrobeCensus failure is deliberately **non-fatal**:
> if it fails for a sample (reads shorter than 50 bp, too few marker-gene
> hits in very small datasets, or an estimate rejected by the module's
> 0.5–20 Mb plausibility check), the run continues and that sample's `cpge`
> fields stay empty (TPM-only fallback). `summary/ags_and_ge.tsv` records the
> per-sample estimates, with `status = unavailable` marking exactly the
> samples that fell back. Note MicrobeCensus estimates need a few hundred
> thousand reads to be accurate — on very small datasets the value is
> mechanical, not meaningful.

> **Read-level tables (`read_*` files, `--read_level_functional`):** the read
> branch uses the same long schema as `function_abundance.tsv` with
> `source = reads` and `backend` naming the tool, but is written to its own
> files and never merged with the contig branch. `abundance_native` is the
> backend's own unit — for Woltka, reads assigned (mates counted separately;
> fractional when a read is divided 1/k across its hit ORFs) with
> `native_unit = reads`. `abundance_tpm` is empty (TPM is contig-branch only).
> `abundance_cpge` uses the same definition as above, with RPK computed over the
> WoLr2 reference ORF lengths: `CPGE = (reads / (ORF length / 1000)) /
> genome_equivalents`, blank when MicrobeCensus is unavailable for the sample
> (`read_sample_summary.tsv` then shows `cpge_status = unavailable`). Read-branch
> and contig-branch values are **not directly comparable**: they count different
> things (reads hitting reference genomes vs reads mapped back to the sample's own
> assembled genes), and read-level profiling recovers the unassembled fraction
> but over-predicts, while the assembly-based branch is more precise. Woltka COG
> accessions are ortholog
> groups (`COG0604`), not the single-letter categories of the contig branch;
> Pfam accessions keep their version suffix as in WoLr2 (`PF00004.32`); MetaCyc
> accessions are pathways (`PWY-5101`); EC numbers are derived via KO.
>
> For SUPER-FOCUS, `abundance_native` is reads with an accepted SEED hit (mates
> counted separately; each read contributes 1 in total, fractional when divided
> 1/k across its best-hit assignments), summed per level into `seed_level1`,
> `seed_level2` and `seed_level3` with path-qualified accessions (`L1`,
> `L1 | L2`, `L1 | L2 | L3`) and the level's own name as `description`; every
> level sums to the reads with a hit. `abundance_cpge` is always empty (a SEED
> hit has no gene length, so no RPK) and `cpge_status = not_applicable`; only
> the three `_native` wide matrices are written. The raw SUPER-FOCUS tables are
> kept for provenance (tab-separated despite the `.xls` suffix), but their `%`
> columns and per-level files are **not** used: SUPER-FOCUS's own level files
> aggregate by level name, merging the level-2 placeholder `-` across 33
> level-1 categories; the pipeline recomputes the levels from the
> function-level counts.

### Database Storage Directory

By default, databases are stored separately from results at `<output_dir>/../databases/`:

```
databases/                            # Database storage (configurable via --databases_dir)
├── bowtie_index/                     # Combined Bowtie2 decontamination index (host + phiX)
├── kraken/                           # Kraken2 database (selected via --kraken2_db)
├── sourmash/                         # Sourmash database
├── taxdump/                          # NCBI taxonomy dump
├── blast/                            # NCBI NT BLAST database
├── deeparg_db/                       # DeepARG database
├── rgi/                              # CARD database for RGI
├── checkm2/                          # CheckM2 database
├── gtdbtk/                           # GTDB-TK database
├── eggnog/                           # eggNOG 7 data for eggNOG-mapper v3 (~44 GB)
├── dbcan/                            # dbCAN database for run_dbcan v5 (~7.4 GB)
├── cog/                              # NCBI COG definitions table (cog-24.def.tab, ~410 KB)
├── bakta/                            # Bakta database v6.0 (full 31.9 GB / light 1.3 GB download)
├── woltka/                           # Web of Life WoLr2 subset for Woltka (~94 GB)
└── superfocus/                       # SUPER-FOCUS DB_90 root for the selected aligner (~1.9-2.5 GB unpacked)
```

The KARGA and KARGVA reference FASTAs are small and staged directly into the work
directory when needed; they are not stored under `databases/`.

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
| `07_functional_annotation/gene_calling/` | `contig_tax_and_arg=true` or `contig_level_functional=true` | Pyrodigal ORF predictions on contigs |
| `07_functional_annotation/eggnog/` | `contig_level_functional=true` | eggNOG-mapper functional annotation of predicted proteins |
| `07_functional_annotation/dbcan/` | `contig_level_functional=true` and `functional_cazy=true` | run_dbcan CAZy annotation with per-tool calls and substrate predictions |
| `07_functional_annotation/gene_abundance/` | `contig_level_functional=true` | featureCounts per-gene read counts (per sample in both assembly modes) |
| `07_functional_annotation/microbecensus/` | `microbecensus=true` and (`contig_level_functional=true` or `read_level_functional != 'none'`) | MicrobeCensus average genome size and genome equivalents per sample |
| `07_functional_annotation/summary/` | `contig_level_functional=true` | Study-level tables: gene/function abundance with TPM and CPGE, wide matrices per ontology, annotated fraction, AGS summary |
| `07_functional_annotation/mags/` | `mag_level_functional=true` (needs `include_binning=true` and ≥2 binners) | Bakta per-bin MAG annotation (GFF3, GBFF, FAA, FNA, TSV, summary per refined bin) |
| `07_functional_annotation/reads/woltka/` | `read_level_functional='woltka'` | Woltka read-level ORF and function tables per sample (no assembly needed) |
| `07_functional_annotation/reads/superfocus/` | `read_level_functional='superfocus'` | SUPER-FOCUS SEED level 1-3 tables and raw SUPER-FOCUS outputs per sample (no assembly needed) |
| `07_functional_annotation/summary/read_*.tsv` | `read_level_functional != 'none'` | Study-level read-branch tables (same schema, `source=reads`), reported separately from the contig branch |
| `07_functional_annotation/contigs/` | `contig_level_metacerberus=true` | MetaCerberus functional annotation results |

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
    --custom_decontamination_index /shared/db/bowtie_index \
    --custom_kraken_db /shared/db/kraken2_standard8 \
    --custom_checkm2_db /shared/db/checkm2/uniref100.KO.1.dmnd \
    --custom_gtdbtk_db /shared/db/release232 \
    -profile singularity
```

### Example 9: Non-Human Host Analysis

For mouse gut microbiome:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --custom_host_fasta /path/to/mouse_genome.fa \
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
    --custom_gtdbtk_db /shared/databases/bugbuster/release232 \
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
    custom_gtdbtk_db              = '/shared/db/release232'
    custom_decontamination_index  = '/shared/db/bowtie_index'
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

The dedicated [`troubleshooting.md`](troubleshooting.md) guide covers
installation, input, resource, container, database, decontamination, cloud and
output issues in depth. The most common cases:

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
# Pre-pull containers (look up the exact tag in the module's main.nf, e.g. modules/nf-core/fastp/main.nf)
singularity pull docker://quay.io/biocontainers/fastp:<tag>

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

# Force fresh run: simply omit -resume
nextflow run main.nf ...

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

*BugBuster v1.1.0dev - Built with Nextflow*
