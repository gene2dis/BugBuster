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
- **Taxonomic Profiling**: Species identification using Kraken2, Sourmash or Sylph
- **Genome Assembly**: Per-sample or co-assembly using MEGAHIT, with optional BBNorm normalization of the assembler input
- **Metagenomic Binning**: MAG recovery with SemiBin (default), MetaBAT2 and COMEBin
- **Bin Refinement**: Quality improvement with Binette, scored with CheckM2 (when two or more binners run)
- **Bin Quality Assessment**: Completeness and contamination with CheckM2
- **Taxonomic Classification**: Bin taxonomy with GTDB-TK
- **ARG Prediction**: At read, contig, and bin levels
- **Functional Annotation**: Contig-level with eggNOG-mapper v3 (Pyrodigal gene calling + eggNOG 7 orthology transfer), run_dbcan v5 CAZy annotation with substrate prediction, per-gene abundance quantification (featureCounts), average genome size estimation (MicrobeCensus) and study-level TPM and copies-per-genome-equivalent tables per functional ontology (KO, COG, EC, Pfam, CAZy; see the container note under the feature toggles); MAG-level with Bakta (per-bin annotation of Binette-refined bins); optional read-level profiling with Woltka against the Web of Life (WoLr2) genomes (KO, EC, COG, Pfam, MetaCyc pathway read counts and copies per genome equivalent) or with SUPER-FOCUS against its SEED subsystem database (SEED subsystem levels 1–3 read counts), or with HUMAnN 4 (alpha 4.0.0a2, MetaPhlAn 4.1.2 prescreen; MetaCyc pathway, KO and EC abundances in RPK with copies per genome equivalent), one backend per run, reported separately from the contig branch; needs no assembly

### Pipeline Workflow

```
Reads → QC (FastP) → Host Removal (bwa-mem3) → Taxonomy (Kraken2/Sourmash/Sylph)
                                            ├──→ Read-level functions (Woltka/WoLr2, SUPER-FOCUS/SEED or HUMAnN 4, optional)
                                            ↓
                                      Assembly (MEGAHIT)
                                            ↓
                                      Contig Filtering (BBMap)
                                            ↓
                              ┌─────────────┴─────────────┐
                              ↓                           ↓
                      Binning (SemiBin;           Contig Analysis
                       + MetaBAT2/COMEBin)         (Metabuli, BlobTools, DeepARG;
                                                    geNomad, optional)
                              ↓
                      Refinement (Binette, ≥2 binners)
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
| QC + Taxonomy only | 16 GB (75 GB to build the human bwa-mem3 index, see below) | 4 | 50 GB |
| QC + Taxonomy + Assembly | 64 GB | 16 | 100 GB |
| Full pipeline with binning | 128 GB | 16 | 200 GB |
| Large datasets / co-assembly | 256 GB | 32 | 500 GB |

- **CPU:** the default aligner, bwa-mem3, needs an x86-64 CPU with AVX2
  (x86-64-v3, Haswell 2013 or newer). On older CPUs use `--read_aligner bowtie2`.
- **Set `--max_memory` to the node's RAM** (default 128 GB). Task requests are
  capped by it, and on a local machine it is also the memory pool Nextflow
  schedules from, so it decides what runs, not the table above. Decontamination
  and per-sample MEGAHIT request 72 GB each on the first attempt
  (`process_high`), although decontamination peaked at ~16 GB. Co-assembly
  MEGAHIT requests memory from its input: 6 GB per GB of gzipped pooled reads,
  at least 72 GB (a co-assembly of 404M clean reads needed ~86 GB at k=21 for
  the graph build alone; for its ~25 GB of reads, estimated, the request is
  ~150 GB). The request is capped by `--max_memory`, so **a large co-assembly
  needs `--max_memory` above its 128 GB default**; otherwise it gets 128 GB
  and may run out of memory (see
  [Out of memory](troubleshooting.md#out-of-memory)). Binette bin refinement
  requests 36 GB on the first attempt (peaks of 11-13 GB measured).
- **The bwa-mem3 index build** needs ~75 GB for the human host (the request is
  ~90 GB for the default gzipped genome, capped by `--max_memory`). It runs
  once per `--databases_dir`: later runs find the stored index
  (`<databases_dir>/bwamem3_index/phiX174_human/contaminants_index`) and skip
  the build. An index built from `--custom_host_fasta` / `--custom_phiX_fasta`
  is rebuilt every run; pass its copy
  (`<databases_dir>/bwamem3_index/contaminants_index`) with
  `--custom_decontamination_index`. Below ~75 GB, pass a pre-built index or use
  `--read_aligner bowtie2` (~6 GB to build).

#### Functional annotation branches

The functional branches add a few memory-heavy steps. Peak memory (RSS) and run
time observed on the pipeline's own validation runs:

| Step | Branch | Peak memory | Run time | Data |
|------|--------|-------------|----------|------|
| `EGGNOG_MAPPER_SEARCH` | contig | **79.9 GB** | 2 h 50 min | real 5-sample co-assembly, 67,812 proteins |
| `EGGNOG_MAPPER_SEARCH` | contig | 23–28 GB | 48–60 min | 2-sample test data (loading the eggNOG DIAMOND database dominates, so even tiny inputs take this long) |
| `EGGNOG_MAPPER_ANNOTATE` | contig | 35.1 GB | 18 min | real co-assembly |
| `RUN_DBCAN` | contig | 5.6 GB | 40 min | real co-assembly |
| `WOLTKA_ALIGN` | read (woltka) | 64.8 GB | — | test data; the WoLr2 index alone needs ≥ 68 GB RAM |
| `HUMANN` | read (humann) | 23.5 GB | ~9 min per sample | test data, full ChocoPhlAn |
| `SUPERFOCUS` | read (superfocus) | ~16 GB | 66–93 min per sample | real samples, DIAMOND |
| `BAKTA_BAKTA` | MAG | 7.9–12.7 GB | 16–24 min per bin | real bins |

These are observations, not limits: needs grow with input size and community
complexity. The eggNOG search requests 16 CPUs and 96 GB on the first attempt;
an out-of-memory failure is retried with 2 × 96 GB, capped at `--max_memory`
(default 128 GB). The 79.9 GB peak above was measured under apptainer, which
does not enforce task memory on a local machine, so make at least ~96 GB
available to it where memory is enforced (docker, SLURM, cloud batch). The annotate step (`process_medium`, 36 GB on
the first attempt) peaked at 35.1 GB, close to that request, so on larger
co-assemblies it may also need its retry.

#### Contig taxonomy (Metabuli) and geNomad

`--contig_tax_and_arg` classifies the contigs with Metabuli against
`refseq_standard` (`METABULI_CLASSIFY`, `process_medium`): in the Phase 9c
benchmark it used 5-12 GB with `--max-ram 32` (the module passes 90 % of the
task's memory) and 443-942 CPU-s per assembly at 16 threads. The database is
110 GB on disk (an 80 GB download, streamed into place: no extra disk for the tarball).
`--genomad_prediction` (geNomad 1.13.0, `process_medium`) needed 17.4-18.4 GB
and 1.1-5.5 CPU-h per assembly. Contig taxonomy no longer uses BLAST against nt
(201,739 CPU-s for one 36.6 Mbp assembly in Phase 9a); see
`docs/benchmarks/contig_taxonomy_2026.md` and
`docs/benchmarks/contig_classifiers_2026.md`.

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

BugBuster automatically downloads required databases on first use and stores them in `<output_dir>/../databases/` by default (configurable via `--databases_dir`), separate from the results directory. Each one goes into a folder named after the database and its registry key (`<databases_dir>/<database>/<key>/`, e.g. `kraken/standard-8/k2_ref_db`; the full layout is under [Database Storage Directory](#database-storage-directory)). A later run with the same `--databases_dir` finds them there and skips the download or build (except the sourmash k-mer database, below) (Nextflow `storeDir`; the log shows `[skipping] Stored process`), whatever its work directory. Selecting another version (e.g. `--kraken2_db gtdb_220`) downloads it into its own folder. The files are moved out of the work directory, so each database is stored once.

- **Refreshing a database**: delete its `<databases_dir>/<database>/<key>/` folder; the next run downloads it again. Do the same if a download was interrupted while being stored (a folder with missing files is taken as complete). RGI's default `--rgi_card_version latest` is stored as `rgi/latest_wildcard_k61`: the first run fixes which CARD release `latest` means for that `--databases_dir` (the stored `versions.yml` records the real version); delete the folder or pin `--rgi_card_version` to update it.
- **Builds from your own files** (the decontamination index from `--custom_host_fasta` / `--custom_phiX_fasta`, RGI with `--custom_rgi_wildcard`, the sourmash taxonomy from a `--custom_sourmash_db` lineages file) have no version to key a folder on: they run in every new run and are copied to `<databases_dir>` for reuse through the matching `--custom_*` option.
- **`-stub` runs** store their placeholder files under `<databases_dir>/stub/`, never where a real run would take them for a database.
- **`--custom_*` paths** are used in place and never copied.
- **The sourmash k-mer database** (`gtdb-reps-rs220-k31.zip`, `--taxonomic_profiler sourmash`) is not downloaded by a pipeline task but staged by Nextflow from its URL into each work directory: a new work directory fetches it again. Download it once and pass it with `--custom_sourmash_db` to avoid that.
- **One run at a time fills a new database folder.** Nextflow has no lock on these folders: two runs started together against an empty `--databases_dir` both download, and a run started while another is still moving a large database into place can take the half-moved folder as complete. Let the first run store its databases before starting others that share the directory.
- **Earlier releases** published databases without the `<key>/` level (e.g. `kraken/k2_ref_db`); this release does not find them there. Pass them with the `--custom_*` options, or delete them and let the pipeline download again.

### Automatic Download Databases

| Database | Size | Used For | Download Trigger |
|----------|------|----------|------------------|
| **phiX174 genome (FASTA)** | 5.4 kB | PhiX contamination removal | `quality_control=true` |
| **Human host genome T2T-CHM13v2.0 (FASTA)** | 940 MB | Host read removal | `quality_control=true` |
| **Kraken2 Standard-8** | 7.5 GB | Taxonomic profiling | `taxonomic_profiler='kraken2'` |
| **Kraken2 GTDB r220** | 497 GB | Taxonomic profiling | `kraken2_db='gtdb_220'` |
| **Sourmash GTDB r220** | 17 GB | Taxonomic profiling | `taxonomic_profiler='sourmash'` |
| **Sylph GTDB r232 (c200 .syl2db) + sylph-tax metadata** | 19 GB + 5.5 MB | Taxonomic profiling | `taxonomic_profiler='sylph'` |
| **Metabuli refseq_standard** | 80 GB download, 110 GB unpacked | Contig taxonomy; its own NCBI taxdump (2026-04-03) feeds BlobTools | `contig_tax_and_arg=true` |
| **geNomad 1.9** | 0.84 GB download, 2.2 GB unpacked | Plasmid / virus / provirus classification of contigs | `genomad_prediction=true` |
| **DeepARG** | 4.8 GB | Contig and bin ARG prediction | `contig_tax_and_arg=true` or `arg_bin_clustering=true` |
| **CARD (RGI)** | 500 MB - 50 GB | AMR gene prediction with pathogen-of-origin | `rgi_prediction=true` |
| **CARD 4.0.1 (K-MARVEL)** | 4.6 MB download, 4.1 MB stored | Read-level ARG assembly with K-MARVEL | `kmarvel_prediction=true` |
| **CheckM2** | 2.9 GB | Bin quality assessment | `include_binning=true` |
| **GTDB-TK r232** | ~61 GB download | Bin taxonomic classification (GTDB-Tk 2.8.0; only R232 data works) | `include_binning=true` |
| **eggNOG 7 (emapper-3.0)** | 44 GB | Contig functional annotation | `contig_level_functional=true` |
| **dbCAN (db_v5-2-9_5-5-2026)** | 7.4 GB | CAZy annotation of predicted proteins | `contig_level_functional=true` (unless `functional_cazy=false`) |
| **NCBI COG 2024 definitions (cog-24.def.tab)** | 410 KB | Maps eggNOG's and WoLr2's COG ids to COG functional categories | `contig_level_functional=true` or `read_level_functional=woltka` |
| **Bakta DB v6.0 full** | 31.9 GB download | MAG (bin) annotation | `mag_level_functional=true` |
| **Bakta DB v6.0 light** | 1.3 GB download | MAG (bin) annotation, reduced annotation sources | `mag_level_functional=true` with `bakta_db='v6.0-light'` (explicit choice, recorded in provenance) |
| **Web of Life WoLr2** | ~94 GB (Bowtie2 index 93.6 GB + coordinates/maps ~0.6 GB) | Read-level functional profiling (Woltka); alignment needs ≥ 68 GB RAM | `read_level_functional='woltka'` |
| **SUPER-FOCUS DB_90** | ~0.74 GB download (~1.9 GB unpacked) for DIAMOND, or ~0.9 GB (~2.5 GB unpacked) for MMseqs2; only the selected aligner's archive | Read-level SEED subsystem profiling (SUPER-FOCUS) | `read_level_functional='superfocus'` |
| **HUMAnN 4 databases + MetaPhlAn vOct22** | ~71 GB download with the full ChocoPhlAn (44.8 GB), ~33 GB with the EC-filtered one (6.9 GB); plus UniRef90 EC-filtered 0.94 GB, utility mapping 2.8 GB and MetaPhlAn `mpa_vOct22_CHOCOPhlAnSGB_202403` with its Bowtie2 index ~22.5 GB | Read-level profiling with HUMAnN 4.0.0a2 | `read_level_functional='humann'` (`humann_db='v4_alpha-full'` or `'v4_alpha-ec_filtered'`) |

### Manual Database Download

For shared HPC systems or repeated runs, pre-download databases to shared storage:

```bash
# Create database directory
mkdir -p /shared/databases/bugbuster/kraken2_standard8

# Download Kraken2 Standard-8
wget -O - https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz | \
    tar -xzf - -C /shared/databases/bugbuster/kraken2_standard8

# Download Sourmash GTDB
wget -P /shared/databases/bugbuster/sourmash/ \
    https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs220/gtdb-reps-rs220-k31.zip
wget -P /shared/databases/bugbuster/sourmash/ \
    https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs220/gtdb-rs220.lineages.reps.csv

# Download Sylph GTDB r232 (19 GB) + its sylph-tax taxonomy metadata (5.5 MB)
# -> --custom_sylph_db .../gtdb-r232-c200-dbv2.syl2db
#    --custom_sylph_taxonomy .../gtdb_r232_metadata.tsv.gz
wget -P /shared/databases/bugbuster/sylph/ \
    https://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r232-c200-dbv2.syl2db
wget -P /shared/databases/bugbuster/sylph/ \
    https://zenodo.org/records/19646381/files/gtdb_r232_metadata.tsv.gz

# Download CheckM2
wget -O /shared/databases/bugbuster/checkm2_db.tar.gz \
    "https://zenodo.org/records/14897628/files/checkm2_database.tar.gz?download=1"
# The tarball holds CheckM2_database/uniref100.KO.1.dmnd; move the .dmnd up so
# it is .../checkm2/uniref100.KO.1.dmnd (--custom_checkm2_db)
mkdir -p /shared/databases/bugbuster/checkm2
tar -xzf /shared/databases/bugbuster/checkm2_db.tar.gz -C /shared/databases/bugbuster/checkm2
find /shared/databases/bugbuster/checkm2 -mindepth 2 -name '*.dmnd' -exec mv {} /shared/databases/bugbuster/checkm2/ \;

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

# HUMAnN 4.0.0a2 database ROOT for the HUMAnN read-level backend (~71 GB with
# the full ChocoPhlAn; for the EC-filtered one use chocophlan_EC_FILTERED.v4_alpha.tar.gz).
# Pass the root as --custom_humann_db. The MetaPhlAn database MUST be
# mpa_vOct22_CHOCOPhlAnSGB_202403 (HUMAnN 4.0.0a2 refuses any other).
H=https://huttenhower.sph.harvard.edu/humann_data
M=https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases
R=/shared/databases/bugbuster/humann
mkdir -p $R/chocophlan $R/uniref $R/utility_mapping $R/metaphlan
wget -O - $H/chocophlan/chocophlan.v4_alpha.tar.gz | tar -xz -C $R/chocophlan
wget -O - $H/uniprot/uniref_ec_filtered/uniref90_annotated_v4_alpha_ec_filtered.tar.gz | tar -xz -C $R/uniref
wget -O - $H/full_mapping_v4_alpha.tar.gz | tar -xz -C $R/utility_mapping
for t in mpa_vOct22_CHOCOPhlAnSGB_202403 bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202403_bt2; do
    wget -O $(basename $t).tar $M/$t.tar && wget -O - $M/$t.md5 | md5sum -c - \
        && tar -xf $(basename $t).tar -C $R/metaphlan; done
# the .pkl and all six .bt2l files must sit directly in $R/metaphlan
find $R/metaphlan -mindepth 2 -type f -exec mv -t $R/metaphlan {} +

# NCBI COG 2024 definitions table (~410 KB; maps eggNOG's COG ids to COG
# functional categories in the contig branch)
wget -P /shared/databases/bugbuster/cog/ \
    https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.def.tab

# Download human host genome (T2T-CHM13v2.0); the pipeline builds the
# combined phiX + host index (bwa-mem3, or Bowtie2 under --read_aligner bowtie2)
# from FASTA on first use
wget -P /shared/databases/bugbuster/ \
    https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
```

### Using Custom Databases

Specify custom database paths to skip automatic downloads:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --custom_sylph_db /shared/databases/bugbuster/sylph/gtdb-r232-c200-dbv2.syl2db \
    --custom_sylph_taxonomy /shared/databases/bugbuster/sylph/gtdb_r232_metadata.tsv.gz \
    --custom_host_fasta /shared/databases/bugbuster/chm13v2.0.fa.gz \
    --custom_checkm2_db /shared/databases/bugbuster/checkm2/uniref100.KO.1.dmnd \
    --custom_gtdbtk_db /shared/databases/bugbuster/release232 \
    --custom_bakta_db /shared/databases/bugbuster/bakta/db \
    --custom_woltka_db /shared/databases/bugbuster/wol2 \
    --custom_superfocus_db /shared/databases/bugbuster/superfocus \
    --custom_humann_db /shared/databases/bugbuster/humann \
    --custom_cog_db /shared/databases/bugbuster/cog/cog-24.def.tab \
    -profile docker
```

### Using Custom Host Genomes

To filter reads from non-human hosts, point the pipeline at the host genome
FASTA (plain or gzipped) — it builds the combined phiX + host index for
`--read_aligner` automatically:

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

If you already have a pre-built combined (phiX + host) index, pass its
directory with `--custom_decontamination_index` instead: a bwa-mem3 or bwa-mem2
index (`*.bwt.2bit.64`) under the default aligner, or a Bowtie2 index (`*.bt2`)
under `--read_aligner bowtie2`.

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
uses the paired reads only (`rgi bwt` accepts a single read pair); Sylph
profiles the paired reads only too (sylph would treat a singleton file as a
separate sample). Sourmash, K-MARVEL (`--kmarvel_prediction`, which pools R1,
R2 and singletons into one read set, as upstream asks for paired-end data) and
assembly use singletons as well.

### Input Requirements

- **Format**: FASTQ or compressed FASTQ (.gz)
- **Paired-end**: Required (R1 and R2)
- **Naming**: Sample names must be unique
- **Paths**: Use absolute paths for cloud/HPC execution
- **Read file names in co-assembly without QC** (`--assembly_mode coassembly
  --quality_control false`): the pooled reads are sorted into R1, R2 and
  singleton lists by file name (`bin/read_role.sh`), so each name needs a role.
  The name minus its extension is split on `_` and `.`; the last `R1`, `R2`,
  `singleton`, `singletons` or `unpaired` token wins (any case), otherwise the
  last bare `1` or `2` (`sample_1.fastq.gz`). A trailing chunk number does not
  override an explicit token (`S1_R1_2.fastq.gz` is R1). A file whose name has
  no role (`sample.fastq.gz`) stops the run with an error. With QC on the clean
  reads are renamed `<sample>_{R1,R2,Singleton}_clean.fastq.gz`, so any input
  names work

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
| `--publish_dir_mode` | `copy` | How to save results: `copy`, `copyNoFollow`, `symlink`, `rellink`, `link`, `move` (see below) |

**Saving space with `--publish_dir_mode link`.** By default every result is copied
from the work directory into `--output`, so large outputs (clean reads with
`--store_clean_reads`, RGI BAMs, bins) take their space twice. With `link`, each
result is a hard link to the file in the work directory: publishing is instant and
takes no extra space. A hard link is a full file, not a pointer:

- copying the results folder elsewhere (`cp -r`, `rsync`, `scp`) copies full files;
- deleting the work directory (by hand, `nextflow clean` or the `low_disk` profile)
  keeps the results; the space is freed when the last of the two names is gone;
- do not edit a published file in place: it is the same file as the one in the work
  directory, which `-resume` reuses.

Hard links only work within one filesystem. If the work directory and `--output` are
on different filesystems (`df <work dir> <output dir>` shows two devices), the first
task that publishes fails with `Failed to publish file ... [link]` (`Invalid
cross-device link` in `.nextflow.log`) and the run stops: Nextflow does not fall back
to copying (tested with Nextflow 25.10.7 and 26.04.4). Use `copy` there, or move the
work directory with `-work-dir`; an interrupted run continues with `-resume`. This is
why `copy` stays the default. `symlink` and `rellink` are different: they publish
links that point into the work directory, which break when it is deleted, and a copy
of the results folder copies only the links unless the copy follows them. Don't use
`move`: it removes outputs that later tasks still read from the work directory and
breaks `-resume`.

### 6.2 Pipeline Execution Options

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `--quality_control` | `true` | `true`, `false` | Enable QC and host filtering |
| `--assembly_mode` | `assembly` | `assembly`, `coassembly`, `none` | Genome assembly strategy |
| `--taxonomic_profiler` | `sylph` | `kraken2`, `sourmash`, `sylph`, `none` | Taxonomic profiling tool (default `sourmash` before 2.0.0) |
| `--include_binning` | `false` | `true`, `false` | Enable binning and refinement |
| `--binners` | `semibin` | `comebin`, `semibin`, `metabat2` | Comma-separated binners to run; ≥2 enables Binette refinement |
| `--args_oap_prediction` | `false` | `true`, `false` | Read-level ARG prediction with ARGs-OAP against SARG (see 6.12; `--read_arg_prediction` is a deprecated alias) |
| `--rgi_prediction` | `false` | `true`, `false` | AMR prediction with pathogen-of-origin (RGI/CARD) |
| `--kmarvel_prediction` | `false` | `true`, `false` | Read-level ARG assembly against CARD protein models with K-MARVEL (needs no assembly; see 6.11) |
| `--contig_tax_and_arg` | `false` | `true`, `false` | Contig taxonomy (Metabuli; BlobTools tables and plots) and ARG prediction (DeepARG); see 6.9 |
| `--genomad_prediction` | `false` | `true`, `false` | Plasmid / virus / provirus classification of the contigs with geNomad (needs an assembly); see 6.9 |
| `--contig_level_functional` | `false` | `true`, `false` | Contig functional annotation (Pyrodigal + eggNOG-mapper v3 + run_dbcan CAZy + featureCounts gene abundance + TPM/CPGE summary tables); see engine note below |
| `--microbecensus` | `true` | `true`, `false` | MicrobeCensus average genome size for CPGE normalization (with the contig or read functional branch; failure is non-fatal — affected samples get no CPGE: TPM-only in the contig tables, native read counts only in the read tables) |
| `--functional_cazy` | `true` | `true`, `false` | run_dbcan v5 CAZy annotation of the predicted proteins (only with the functional branch; reported separately from the eggNOG CAZy calls, never merged) |
| `--dbcan_consensus` | `recommended` | `recommended`, `any` | Which dbCAN calls feed the summary tables (`recommended` = supported by ≥2 tools; `any` = per-tool union) |
| `--mag_level_functional` | `false` | `true`, `false` | MAG-level functional annotation: Bakta on every refined bin (requires `--include_binning` and ≥2 `--binners`; see the note below) |
| `--read_level_functional` | `none` | `woltka`, `superfocus`, `humann`, `none` | Read-level functional profiling backend (needs no assembly; see the notes below) |
| `--woltka_uniq` | `false` | `true`, `false` | Woltka: leave reads whose reported alignments hit several ORFs unassigned instead of dividing them 1/k |
| `--superfocus_aligner` | `diamond` | `diamond`, `mmseqs2` | SUPER-FOCUS search backend (DIAMOND blastx or MMseqs2); selects which DB_90 archive is downloaded |
| `--arg_bin_clustering` | `false` | `true`, `false` | ARG clustering for HGT inference: the assembly's DeepARG hits split by refined bin (no per-bin DeepARG run), then MMseqs2 |
| `--min_read_sample` | `0` | Integer ≥ 0 | Minimum reads required after QC |
| `--store_clean_reads` | `false` | `true`, `false` | Publish the decontaminated reads to `<output>/clean_reads/` (publishing only; `-resume` is the task cache) |
| `--featurecounts_multimap` | `primary` | `primary`, `all`, `none` | featureCounts multi-mapping policy for gene counts: primary alignments only, every reported alignment (`-M`), or no multi-mapping reads |

> **Container note for `--contig_level_functional`:** the branch uses eggNOG-mapper
> 3.0.0-beta6, whose authors publish only an Apptainer image. Under
> `-profile singularity` or `-profile apptainer` the pipeline uses that official image;
> under docker/podman and the docker-based cloud profiles (`aws`, `gcp`, `azure`) it
> uses `ghcr.io/gene2dis/bugbuster-eggnog-mapper`, an image of the same beta6 built by
> this pipeline from the authors' own recipe, with the same tool versions (both were
> checked to give identical results). Both engines apply a run-time fix for an
> eggNOG-mapper beta6 bug that cost some genes their COG category (upstream issue
> #620; see troubleshooting). The pipeline will switch to the official bioconda image
> once eggNOG-mapper v3.0.0 final is released.

> **Note for `--mag_level_functional`:** Bakta annotates the refined bins, one task
> per bin, publishing per-bin GFF3/GBFF/FAA/FNA/TSV/summary files under
> `07_functional_annotation/mags/<sample>/`. It requires `--include_binning` **and at
> least two `--binners`**: Binette refinement and its completeness/contamination
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
> summarized to KO, EC (via KO), COG functional categories (via KO), Pfam families
> and MetaCyc pathways from the WoLr2 maps, using the contig branch's COG and Pfam
> vocabularies (details below). **A read contributes its full weight once to each distinct
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

> **Note for `--read_level_functional humann`:** HUMAnN 4 is an **alpha** release;
> the pipeline pins exactly HUMAnN 4.0.0a2 with MetaPhlAn 4.1.2 and the MetaPhlAn
> database `mpa_vOct22_CHOCOPhlAnSGB_202403` (the only combination this HUMAnN
> build accepts — an existing MetaPhlAn 4 database such as vJun23 is refused).
> Each sample's host-removed R1, R2 and singleton reads are concatenated into one
> input (mates counted as separate reads). HUMAnN prescreens the community with
> MetaPhlAn, aligns the reads to the pangenomes of the detected species
> (ChocoPhlAn, nucleotide search), then searches the remaining reads against
> UniRef90 (translated search; HUMAnN 4 offers only the EC-filtered UniRef90
> database). It runs with `--count-normalization RPKs`, so **`abundance_native` is
> RPK** (`native_unit = rpk`) and copies per genome equivalent are computed as for
> Woltka (`CPGE = RPK / genome_equivalents`, blank without MicrobeCensus). Three
> ontologies are reported: `metacyc` (HUMAnN's MetaCyc pathway abundance), and
> `ko` / `ec`, which the pipeline derives from HUMAnN's gene families with the
> HUMAnN 4 mapping files — a gene family carrying several terms contributes its
> full RPK to each (the same intentional double counting as the other branches).
> HUMAnN 4 gene families are a mix of UniRef90 and UniClust90 clusters: only
> UniRef90 families can map to a KO, while EC maps cover both, so KO coverage is
> lower than EC coverage. **KO tables from this backend are sparse:** HUMAnN
> 4.0.0a2's utility mapping ships a legacy UniRef90→KO map that does not cover
> most HUMAnN 4 gene families — on the nf-core test reads only 17–18 % of the
> gene-family abundance mapped to a KO, against 39–61 % for EC (UniClust90
> families were only 3–4 % of the abundance, so they are not the cause). Check
> the `ko` row of `read_annotated_fraction.tsv` for your data before relying on
> HUMAnN KO tables; for KO profiles, the Woltka backend or the contig branch are
> better sources. `read_annotated_fraction.tsv` reports, for HUMAnN, the
> share of reads mapped by HUMAnN (`any` row: reads given minus `READS_UNMAPPED`)
> and, per ontology, the share of the gene-family RPK that carries a term (blank
> read columns; for `metacyc`, the share of pathway abundance outside `UNMAPPED` /
> `UNINTEGRATED`). The raw HUMAnN tables (MetaPhlAn profile, gene families,
> MetaCyc reactions, pathway abundance, log) are published per sample; HUMAnN
> 4.0.0a2 writes no pathway-coverage table. Two defects of this alpha are worked
> around: its own `humann_renorm_table` / `humann_regroup_table` mishandle the
> `READS_UNMAPPED` row (the pipeline does its own regrouping — drop that row if you
> post-process the raw tables), and a sample name containing `s__` or `t__`
> crashes its MetaPhlAn-profile parsing, so such names are rejected at launch.
> Database footprint: ~71 GB with the full ChocoPhlAn (`--humann_db
> v4_alpha-full`, default) or ~33 GB with the EC-filtered one
> (`v4_alpha-ec_filtered`: fewer gene families and KOs); the MetaPhlAn prescreen
> loads a ~20 GB Bowtie2 index. The branch runs on any container engine and needs
> no assembly.

> **Combining the functional branches:** `--contig_level_functional`,
> `--mag_level_functional` and `--read_level_functional` are independent and can
> all be enabled in one launch (one read backend per run). MicrobeCensus then
> runs once per sample and its genome equivalents feed both the contig and the
> read tables; the read branch still writes only its own `read_*` files and
> per-backend directories, never the contig-branch tables. With
> `--contig_level_functional` on, all three read backends were verified alongside
> the contig branch.

### 6.3 Database Selection Options

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `--phiX_index` | `phiX174` | `phiX174` | PhiX genome for decontamination |
| `--host_db` | `human` | `human` | Host genome for read removal |
| `--kraken2_db` | `standard-8` | `standard-8`, `gtdb_220` | Kraken2 database version |
| `--sourmash_db` | `gtdb_220_k31` | `gtdb_220_k31` | Sourmash database version |
| `--sylph_db` | `gtdb_r232` | `gtdb_r232`, `gtdb_r226`, `gtdb_r220` | Sylph database (GTDB reps, c200) + matching sylph-tax metadata |
| `--checkm2_db` | `v3` | `v3` | CheckM2 database version |
| `--gtdbtk_db` | `release_232` | `release_232` | GTDB-TK database release (only R232 works with the pinned GTDB-Tk 2.8.0) |
| `--eggnog_db` | `emapper-3.0` | `emapper-3.0` | eggNOG 7 data for eggNOG-mapper v3 |
| `--dbcan_db` | `db_v5-2-9_5-5-2026` | `db_v5-2-9_5-5-2026` | dbCAN database release for run_dbcan v5 |
| `--bakta_db` | `v6.0-full` | `v6.0-full`, `v6.0-light` | Bakta database flavor; the light DB changes annotation results, so selecting it is always explicit and is recorded in provenance (`software_versions.yml`) |
| `--woltka_db` | `wolr2` | `wolr2` | Web of Life release for the Woltka read-level backend |
| `--superfocus_db` | `db90` | `db90` | SUPER-FOCUS database for the SUPER-FOCUS read-level backend (DB_90, figshare CC0) |
| `--humann_db` | `v4_alpha-full` | `v4_alpha-full`, `v4_alpha-ec_filtered` | HUMAnN 4 database set: full (44.8 GB) or EC-filtered (6.9 GB) ChocoPhlAn, plus UniRef90 EC-filtered, utility mapping and MetaPhlAn vOct22 (recorded in provenance) |
| `--cog_db` | `cog-24` | `cog-24` | NCBI COG definitions table mapping eggNOG's (contig branch) and WoLr2's (Woltka read backend) COG ids to COG functional categories |
| `--kmarvel_db` | `card_v4.0.1` | `card_v4.0.1` | CARD release for K-MARVEL (CARD 4.0.1 data tarball, the release upstream's `get_ref_db.sh` fetches) |
| `--metabuli_db` | `refseq_standard` | `refseq_standard`, `gtdb232` | Metabuli database for contig taxonomy (`gtdb232`: GTDB names, 634 GB download, not benchmarked) |
| `--genomad_db` | `v1.9` | `v1.9` | geNomad database release (pinned Zenodo file) |

### 6.4 Custom Database Paths

| Parameter | Description |
|-----------|-------------|
| `--custom_decontamination_index` | Path to a pre-built combined index directory (host + phiX) matching `--read_aligner`: bwa-mem3/bwa-mem2 (`*.bwt.2bit.64`, default) or Bowtie2 (`*.bt2`) |
| `--custom_phiX_fasta` | Path to a custom phiX genome FASTA (index is built by the pipeline) |
| `--custom_host_fasta` | Path to a custom host genome FASTA (index is built by the pipeline) |
| `--custom_kraken_db` | Path to custom Kraken2 database directory |
| `--custom_sourmash_db` | List of paths: `["kmer.zip", "lineages.csv"]` |
| `--custom_sylph_db` | Path to a sylph genome database file (`.syldb` or `.syl2db`); requires `--custom_sylph_taxonomy`; outputs are labelled with its file name |
| `--custom_sylph_taxonomy` | Path to the sylph-tax metadata file of the same release (accession, tab, `d__...;s__...` lineage; `.tsv` or `.tsv.gz`); requires `--custom_sylph_db` |
| `--custom_checkm2_db` | Path to CheckM2 database file |
| `--custom_gtdbtk_db` | Path to the directory directly containing the unarchived GTDB-Tk reference data (e.g. the extracted `release232/`); only R232 data works with the pinned GTDB-Tk 2.8.0, R220/R226 does not (see `docs/parameters.md`) |
| `--custom_deeparg_db` | Path to DeepARG database directory |
| `--custom_metabuli_db` | Path to a Metabuli database directory (`db.parameters`, `diffIdx`, `info`, `taxonomyDB`, and its taxdump `nodes.dmp`/`names.dmp` at most one directory down, e.g. `ncbi-taxdump/`) |
| `--custom_genomad_db` | Path to a geNomad database directory (`genomad_db`, `version.txt`, ...; geNomad 1.13 takes database 1.9) |
| `--custom_rgi_card_db` | Path to pre-prepared CARD database directory |
| `--custom_rgi_wildcard` | Path to WildCARD directory (use with `--custom_rgi_card_db`) |
| `--custom_kmarvel_db` | Path to a K-MARVEL database directory: `CARD_protein_all.fasta`, `metadata.tsv` (CARD `aro_index.tsv`) and `snps.txt` (CARD `snps.txt`); no BLAST index needed. Its `DB_VERSION`, or the FASTA's md5 without one, goes to `software_versions.yml` |
| `--custom_eggnog_db` | Path to eggNOG 7 data directory (emapper-3.0 layout, see `docs/parameters.md`) |
| `--custom_dbcan_db` | Path to dbCAN database directory (run_dbcan v5 layout, see `docs/parameters.md`) |
| `--custom_bakta_db` | Path to Bakta database directory (schema 6 layout, see `docs/parameters.md`) |
| `--custom_woltka_db` | Path to a local WoLr2 mirror in the FTP layout (see `docs/parameters.md`) |
| `--custom_superfocus_db` | Path to a SUPER-FOCUS database ROOT, the directory containing `db/` (see `docs/parameters.md`) |
| `--custom_humann_db` | Path to a HUMAnN database ROOT with `chocophlan/`, `uniref/`, `utility_mapping/` and `metaphlan/` (see `docs/parameters.md`) |
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
| `--sylph_read_seq_id` | unset (automatic) | Sylph read identity (%), e.g. `99.5`; corrects the unknown fraction (Classified, estimated read pairs) |

Sylph runs with its defaults plus `-u`; its abundance tables are described in
`docs/parameters.md` (*Sylph outputs*).

### 6.7 Assembly Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assembly_read_reduction` | `none` | `normalize`: BBNorm on the reads MEGAHIT gets (pooled under co-assembly) |
| `--megahit_k_list` | MEGAHIT's own list | MEGAHIT `--k-list`; `21,41,61,81,101,121,141` for large co-assemblies |
| `--megahit_min_count` | `2` | MEGAHIT `--min-count`; `3` for co-assemblies that run out of memory or time |
| `--bbmap_length` | `1000` | Minimum contig length after BBMap filtering |

**Only the assembler input is reduced.** With `--assembly_read_reduction
normalize`, MEGAHIT assembles the normalized reads, and everything after it
(BBMap's contig filter, the read alignments, coverage, binning and gene
counting) uses the full clean reads, so abundances and bin depths are not
affected. Normalization pays off on deep co-assemblies; on a 19.4 M-pair
co-assembly it cost more CPU time than it saved (see
[the assembly benchmark](benchmarks/assembly_2026.md) and `docs/parameters.md`
for the numbers behind each setting).

### 6.8 Binning Options

#### Binner Selection

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `--binners` | `semibin` | `comebin`, `semibin`, `metabat2` | Comma-separated list of binners to run. At least one required. Binette refinement runs when ≥2 binners are selected. |

**Selection logic:**
- **1 binner selected** → runs that binner only, no refinement, binner output used directly; CheckM2 scores its bins in a separate step
- **≥2 binners selected** → runs selected binners + Binette bin refinement; Binette's own CheckM2 run scores the raw and the refined bins, and the pipeline publishes those reports in `quality/checkm2/` (no separate CheckM2 step)

**Co-assembly coverage.** Under `--assembly_mode coassembly` each sample's reads
are aligned once against the co-assembly. Those per-sample BAMs give MetaBAT2
one depth column per sample, count genes for the functional branch and give
each refined bin its per-sample coverage (`coverage/<sample>_bin_depth.tsv`);
SemiBin2 (whose pre-trained models take one coverage column) and COMEBin get
their `samtools merge`, the pooled BAM. On a five-sample co-assembly the
per-sample depths raised MetaBAT2's MIMAG high-quality bins from 11 to 23
(`docs/benchmarks/binning_2026.md`).

> **COMEBin is costly.** It runs on CPU only (the pipeline has no GPU path; 8 CPUs
> per sample) and took 38 CPU-hours (6 h wall) on a 37 Mbp per-sample assembly,
> about 370 times SemiBin2; its cost grows with the assembly, about 1-1.6
> CPU-hours per Mbp of contigs. In the Phase 7a benchmark
> (`docs/benchmarks/binning_2026.md`) adding it to SemiBin2 + MetaBAT2 gave
> Binette one more genome and one more high-quality bin on one sample, and the
> only passing bin on a small one. A failed attempt is retried with smaller
> batches; if the last attempt fails too, its bin folder holds `FAILED.txt` and
> the run goes on without it.

#### Basic Binning Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metabat_minContig` | `2500` | Minimum contig length for MetaBAT2 |
| `--bin_min_completeness` | `50` | Minimum CheckM2 completeness (%) of a refined bin (used when ≥2 binners; was `--metawrap_completeness`, still accepted with a warning) |
| `--bin_max_contamination` | `10` | Maximum CheckM2 contamination (%) of a refined bin (used when ≥2 binners; was `--metawrap_contamination`, still accepted with a warning) |
| `--semibin_env_model` | `human_gut` | SemiBin environment model |
| `--semibin_coassembly_mode` | `pretrained` | Co-assembly: SemiBin2 pre-trained model on the pooled BAM (`pretrained`) or a model trained on the co-assembly (`train_self`, ~40× the CPU time) |
| `--semibin_coassembly_env_model` | `global` | SemiBin2 pre-trained model for the co-assembly |
| `--semibin_random_seed` | `11` | SemiBin2 `--random-seed`, so reruns give the same bins |
| `--comebin_random_seed` | `11` | COMEBin `-s` seed, so reruns give the same bins |
| `--gtdbtk_pplacer_cpus` | `4` | Threads for GTDB-Tk's pplacer step (capped at the task's CPUs) |

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

### 6.9 Contig Taxonomy, geNomad and DeepARG Options

`Contig_tax_and_arg_prediction.tsv` (`05_arg_prediction/contig_level/summary/`)
has one row per contig, or per ARG gene on a contig with several:

- `contig_length`, `GC`, `coverage`: from BlobTools.
- `superkingdom` ... `species` (with `class`), `taxonomy_rank`, `taxonomy_score`:
  Metabuli's lineage, the rank it called at and its score (`unclassified`, 0
  for no call). `superkingdom` is NCBI's `domain` (renamed in 2025), or
  Viruses. Calls below `--metabuli_min_score` are not made, and species calls
  need `--metabuli_min_sp_score` (else genus or higher).
- `bin`, `bin_gtdb_classification` (with `--include_binning`): the refined bin
  and its GTDB-Tk string, never copied into the NCBI ranks.
- `mobile_class` (`chromosome`, `plasmid`, `virus`), `plasmid_score`,
  `virus_score`, `provirus` (geNomad provirus regions overlapping the ARG gene),
  `virus_taxonomy`, `plasmid_conjugation_genes` (with `--genomad_prediction`).
  geNomad's plasmid class also takes integrative conjugative elements and
  transposon cargo. `00_summary/tables/arg_contig_mobile.tsv` counts the ARG
  genes per class on plasmid and virus contigs and in proviruses.
- `ARG_gene` ... `counts`: DeepARG.

BlobTools' tables and phylum plot use Metabuli's calls and the database's own
taxdump; its 2025 `domain` rank is renamed `superkingdom` for BlobTools 1.1.1,
which otherwise wrote `undef` as every contig's superkingdom.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metabuli_min_score` | `0.2` | Metabuli `--min-score` (0-1) |
| `--metabuli_min_sp_score` | `0.5` | Metabuli `--min-sp-score` (0-1): minimum score for a species call |


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

### 6.11 K-MARVEL Options

Parameters for K-MARVEL read-level ARG detection (`--kmarvel_prediction`).
[K-MARVEL](https://bitbucket.org/amr-avengers/k-marvel) assembles CARD protein
ARGs straight from the reads (translated k-mer search, gene-centric protein
assembly), removes chimeric assemblies, keeps housekeeping-gene variant models
only when they carry a resistance SNP from CARD's `snps.txt`, and confirms
every assembled ARG with blastp. It runs per sample on the host-removed reads
(R1, R2 and singletons pooled) in short-read mode, whose k-mer (6) and
extension (9) sizes are fixed by the tool and not exposed.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--kmarvel_min_coverage` | `90` | Minimum coverage (%) of the CARD reference protein by an assembled ARG |
| `--kmarvel_min_identity` | `90` | Minimum identity (%) of an assembled ARG to its CARD reference protein |
| `--kmarvel_min_read_length` | `100` | Minimum read length (bp); must be below 200 in short-read mode. Reads shorter than this after fastp trimming are ignored |

**Notes**:
- K-MARVEL has no bioconda package; the pipeline runs it from its own image
  (`ghcr.io/gene2dis/bugbuster-kmarvel`, built from `containers/kmarvel/`),
  which also carries BLAST+ 2.17.0.
- The default database is the CARD 4.0.1 data release (~4.6 MB). CARD is free
  for academic use; commercial use needs a licence from CARD.
- A sample with no ARG found still gets a (header-only) table and an all-zero
  column in the summary matrices. If blastp fails for a sample with ARGs, the
  task fails instead of publishing unfiltered results.
- Counts are K-MARVEL's supporting reads, which count read pairs (both mates
  of a fragment count once). They are not normalised here; the results summary
  turns them into RPM and RPKM of clean read pairs.
- The allele K-MARVEL reports among close variants of one gene is not stable
  (simulated TEM-1 reads were called TEM-117, TEM-192 or TEM-246 depending on
  the read subset). Compare samples with `KMARVEL_gene_family_matrix.tsv`;
  treat the per-ARO matrix as allele-level detail.

### 6.12 ARGs-OAP Options

[ARGs-OAP](https://github.com/xinehc/args_oap) 4.0.0 (`--args_oap_prediction`)
profiles ARGs against its structured SARG database. Stage one estimates each
sample's 16S copies (GreenGenes) and cells (essential single-copy marker
genes) and extracts ARG-like reads with DIAMOND; stage two searches those reads
with BLASTX against SARG and normalises the hits by SARG type (antibiotic
class), subtype and gene. Both stages run in one task per sample on the
host-removed reads. SARG ships inside the image, so there is no database to
download or point to.

Outputs: `05_arg_prediction/read_level/args_oap/<sample>_args_oap/` holds
ARGs-OAP's own files (`metadata.txt` with nRead, n16S and nCell,
`normalized_cell.*`, `normalized_16S.*`, `ppm.*` and the other tables by type,
subtype and gene). `args_oap_summary/` holds `ARGS_OAP_REPORT`'s run-level
`args_oap_long.tsv` (one row per sample and SARG type, subtype or gene:
read count, copies, copies per cell, per 16S copy, ppm) and `args_oap_cells.tsv`
(nRead, n16S, nCell, SARG reads and a status per sample).

**Notes**:
- Copies per cell = SARG copies / nCell, where a copy is a read weighted by
  the share of the reference it covers (a read on a two- or three-component
  reference also counts 1/2 or 1/3). It is not RPKM: compare ARGs-OAP with
  K-MARVEL / RGI only qualitatively. ppm counts single reads.
- ARGs-OAP treats a singleton file as a separate sample and normalises it by
  its own cells; `args_oap_summary/` adds each sample's singletons to its pairs
  and recomputes the normalisation over both, so its numbers differ slightly
  from ARGs-OAP's own per-sample tables when a sample has singletons.
- A sample where stage two has nothing to normalise (no 16S or cell estimate:
  stage two is skipped; or no read left after the SARG search) keeps
  `metadata.txt` only, and `stage_two_status.txt` says `not_normalised` or
  `no_sarg_hit` (else `ok`); the report carries that status instead of failing.
  If stage one drops an input file after an error (ARGs-OAP only logs a
  warning), the task fails (exit 137, retried with more memory, when a stage
  one process was killed for memory). An empty singleton file is skipped.
- The image indexes its bundled databases on first use; the task runs a copy
  of the `args_oap` package from its work directory so this also works on a
  read-only singularity/apptainer image (~12 s per task).
- `--read_arg_prediction` is a deprecated alias of `--args_oap_prediction`
  (it used to run KARGA, KARGVA and ARGs-OAP together; KARGA and KARGVA were
  removed).

### 6.13 Read Aligner Options

`--read_aligner` picks the aligner for host/phiX removal and reads-to-contigs
mapping (contig depth, binning, gene counting, co-assembly bin depth).
Woltka and HUMAnN always use bowtie2 with their own databases.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--read_aligner` | `bwamem3` | `bwamem3` (bwa-mem3 0.15.0) or `bowtie2` (fallback) |
| `--decontam_min_score` | `30` | bwa-mem3 minimum alignment score (`-T`) for a contaminant hit, per mate |
| `--bwamem3_min_identity` | `0.95` | Minimum identity of a read-to-contig/bin alignment (`0` = off) |
| `--bwamem3_min_aligned_length` | `50` | Minimum aligned bases of a read-to-contig/bin alignment, or 80% of a shorter read (`0` = off) |
| `--bwamem3_dense_sa` | `false` | Denser suffix array for the decontamination index (faster, ~12 GB more RAM) |

bwa-mem3 replaced bowtie2 as the default after a benchmark on five real
metagenomes ([docs/benchmarks](benchmarks/bwamem3_vs_bowtie2.md)): host removal
4.1× faster, reads-to-contigs mapping 3.9–6.2× faster, the same bins and GTDB
species, and a stricter host removal (a pair is clean only when both mates are
unmapped). Use `--read_aligner bowtie2`:

- on a CPU without AVX2 (x86-64-v3, Haswell 2013+): the `BWAMEM3_*` tasks
  stop with "Illegal instruction" (exit status 132);
- on nodes with less than ~75 GB RAM, which the human bwa-mem3 index build
  needs (bowtie2: ~6 GB, but about an hour). The build runs once; pass the
  published index with `--custom_decontamination_index` in later runs;
- to keep using a Bowtie2 `--custom_decontamination_index` from an earlier
  release. A Bowtie2 index is rejected at launch under `bwamem3`.

Clean-read counts from the two aligners are not directly comparable: bowtie2
keeps a pair without a concordant alignment, so a pair with a single host mate
survives.

#### Bowtie2 options

Advanced parameters for Bowtie2 read alignment during host filtering
(`--read_aligner bowtie2` only):

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

### 6.14 Taxonomy Visualization Options

Control taxonomic output visualization and formatting:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--taxonomy_plot_levels` | `Phylum,Family,Genus,Species` | Comma-separated taxonomic levels to plot |
| `--taxonomy_top_n_taxa` | `10` | Maximum taxa per plot; each sample's most abundant taxa are included first, then the most abundant overall |
| `--create_phyloseq_rds` | `false` | Generate R phyloseq RDS files for downstream analysis |

**Taxonomic Levels:** Domain (D), Phylum (P), Class (C), Order (O), Family (F), Genus (G), Species (S)

### 6.15 Database Storage Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--databases_dir` | `<output>/../databases` | Directory for storing downloaded databases (separate from results) |

By default, databases are stored in a `databases/` directory at the same level as your output directory. Later runs with the same `--databases_dir` reuse them without downloading again; see [Database Management](#4-database-management).

### 6.16 Resource Limit Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `240.h` | Maximum time per process |

### 6.17 Results Summary Options

The run-level summary in `00_summary/` consolidates the QC, taxonomy, read-level ARG (K-MARVEL, RGI, ARGs-OAP),
read functional and MicrobeCensus results: tidy and wide TSV tables, an Excel workbook (`results_summary.xlsx`, with a
README sheet defining every column), PDF/PNG/SVG figures, a self-contained HTML report (`results_summary.html`, no
network access needed to view it), a README (`README.md` and `README.docx`) and `LLM_GUIDE.md`. Each part appears only
when its stage ran; the README and the guide list every pipeline stage and whether it ran in each results directory.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--results_summary` | `true` | Write `00_summary/` |
| `--sample_metadata` | none | TSV/CSV with a `sample` column (or sample names in the first column) and any metadata columns; carried into every table |
| `--summary_group_by` | none | Metadata column that groups samples (subject, site, ...) in the figures and the per-group ARG detection table |
| `--summary_order_by` | none | Metadata column that orders samples within a group (visit, day, ...) in the figures and the per-group detection patterns |
| `--summary_batch_label` | output directory name (its parent's when it is `results`) | Label in the `batch` column; set it when the default does not tell the batches apart |
| `--summary_min_clean_pairs` | `1000000` | Samples with fewer clean read pairs are flagged as low depth; a gene or taxon not found there is reported as `NA_low_depth`, not as absent |
| `--summary_top_n` | `40` | Maximum rows in the summary ARG heatmap |
| `--summary_taxa_per_sample` | `5` | The taxonomy heatmaps show the union of each sample's N most abundant taxa (no row cap), grouped by phylum, A-Z within |
| `--summary_depth_marks` | `10000000,20000000` | Clean read-pair depths marked in the QC figure and counted in the report header |
| `--rgi_summary_min_coverage` | `80` | An RGI gene counts when reads cover at least this % of its best reference sequence (summaries only) |
| `--card_prevalence` | none | WildCARD `card_prevalence.txt` (an RGI WildCARD database has it in `wildcard/`); adds candidate hosts and CARD plasmid association to each ARG call |
| `--arg_priority_table` | `assets/arg_priority_tiers.tsv` | ARG reporting tier rules |
| `--summary_watchlist` | `assets/summary_watchlist.tsv` | Taxa reported in the watchlist table and figure |

How the summary reports the results:

- **Read pairs throughout.** `Reads_report.csv` counts single reads (both mates); the summary divides by two, the unit K-MARVEL and sylph use.
- **ARG abundance is normalised.** K-MARVEL's supporting read pairs become RPM (per million clean read pairs) and RPKM (RPM per kb of the CARD reference gene). The raw matrices in `kmarvel_summary/` are not comparable between samples of different depth.
- **ARG tiers.** Each call is assigned a tier by `--arg_priority_table`: 1 clinical priority (carbapenemases, mcr, vanA/B/D/M, mecA/C, 16S rRNA methyltransferases, cfr/optrA/poxtA, tet(X3+)), 2 other acquired, 3 intrinsic, efflux or regulator, 4 mutation dependent. The rules are regular expressions on CARD gene families and names; copy the file to change them.
- **Mutation-dependent hits are checked.** For variant-model hits (gyrA, parC, rpoB, ...), K-MARVEL lists every difference from the reference. The summary checks them against the CARD `snps.txt` of the K-MARVEL database. Only rows with `resistance_snp_status = card_snp_present` carry a mutation CARD lists for resistance; the others are detections of the gene, not of resistance.
- **Genes, not alleles.** Read-level allele calls among close variants are unstable and share their reads (TEM-1 and TEM-256 with 16,032 and 16,029 pairs), so each gene is represented by its best-supported allele, never by the sum.
- **Candidate hosts are co-occurrence.** With `--card_prevalence`, each ARG call lists the species of the same sample that CARD finds carrying that ARO. Reads are not linked to a host; confirming one needs assembly or long reads.
- **Sylph confidence.** `taxonomy_long.tsv` adds sylph's ANI, ANI interval and coverage per genome; `species_call_uncertain` marks genomes whose ANI interval reaches below 95%.
- **RGI genes need a coverage cutoff.** `rgi bwt` lists every gene with at least one mapped read, fragments of related genes included (548 genes in one batch1 sample). The summary counts an RGI gene when reads cover at least `--rgi_summary_min_coverage` % (80) of its best-covered reference sequence; the raw RGI tables are unchanged. RGI counts mates separately, so its read pairs are mapped reads / 2.
- **RGI pathogen of origin.** RGI's k-mer step assigns the gene's reads to species. `rgi_kmer_hosts` lists the top ones and `rgi_host_check` says whether the top host is a watchlist pathogen. WildCARD holds resistance-gene sequences from many species, gut commensals among them, so a tier 1 gene whose reads point to a commensal (poxtA from *Blautia*, vanD-type genes from *Enterocloster*) is flagged as a possible related gene rather than reported as a pathogen's resistance.
- **Method comparison.** When K-MARVEL and RGI both ran, `arg_method_comparison.tsv`, `arg_method_overlap.tsv`, an UpSet plot and a per-category bar chart show, gene by gene (alleles collapsed), what each method found. A K-MARVEL gene RGI saw below its cutoff is marked `rgi_below_cutoff`, which is partial coverage, not absence. The ARG matrices and figures use one method (`--arg-primary` in the script; K-MARVEL when present); `arg_long.tsv` holds every method's calls.
- **ARGs-OAP (SARG).** With `--args_oap_prediction`, the summary reads `ARGS_OAP_REPORT`'s tables into `arg_sarg_long` (one row per sample and SARG subtype), `arg_sarg_type_cpc` and `arg_sarg_subtype_cpc` (copies per cell, samples as columns), `samples_qc` `args_oap_n16S`, `args_oap_ncell`, `sarg_n_subtypes` and `sarg_copies_per_cell`, and a heatmap by SARG type. SARG names its types and subtypes differently from CARD, so ARGs-OAP is reported in its own section, not in the K-MARVEL vs RGI comparison, and has no tiers. Copies per cell is not RPKM: compare samples within one unit only. A sample ARGs-OAP could not normalise (no 16S or cells) is blank, not 0.
- **Watchlist.** The default `assets/summary_watchlist.tsv` lists the ESKAPE pathogens and the bacteria of the WHO Bacterial Priority Pathogens List 2024 and the CDC Antibiotic Resistance Threats Report 2019, plus a few other species; each row's note names its source. The report's Methods section shows the list.
- **Read functional profile and MicrobeCensus.** With `--read_level_functional`, `function_read_long` holds AGGREGATE_READ_FUNCTIONS' table of every batch and `function_read_<ontology>` one wide table per ontology, in copies per genome equivalent (CPGE) when the backend computes it (Woltka, HUMAnN), else in its native unit; a sample without a MicrobeCensus estimate is `NA_no_ags` there. `samples_qc` adds the reads assigned, the annotated fraction, the MicrobeCensus average genome size and genome equivalents. A COG-category heatmap gives an overview. An ORF can carry several terms, so the tables are not compositional.
- **README and LLM guide.** `README.md` / `README.docx` (written with the Python standard library) list the stages that ran in each results directory, every table and figure with its rows and meaning, the units, the tiers, the caveats and the software versions. `LLM_GUIDE.md` describes the summary folder for a language model or a script: run facts as YAML, the stages that ran, methods, every table column and figure, pitfalls (single reads vs pairs, non-comparable units, NA codes) and recipes. The summary folder is the part meant to be shared, so neither file names a local path: paths are relative to the folder, batches appear by label and input files by name. Which stages ran is read from `software_versions.yml`; a stage that ran without its table in the results directory is reported as "ran, no table found".
- **Attention list and reading guide.** The HTML header lists the tier 1 genes (found by every method, by one only, or with a commensal host) and the low-depth samples. The first section, "Reading this report", explains the units, the ARG tiers (built from the tier rule file in use) and the ARG methods before any result is shown.

The same script summarises several finished results directories (batches) without rerunning any tool; see [Example 13](#example-13-summarise-several-batches).

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
| `low_disk` | Deletes the work dir after a successful run (`cleanup = true`; does not lower peak usage during the run; a completed run cannot be resumed, an interrupted one can) plus `--store_clean_reads` — see [`DISK_OPTIMIZATION.md`](DISK_OPTIMIZATION.md) |
| `test` | Run with minimal test data |

> **`low_disk` with the functional branches:** the profile only cleans the work
> directory. The functional annotation databases are stored at
> `--databases_dir`, outside the work directory (eggNOG 7 ~44 GB, dbCAN ~7.4 GB,
> Bakta up to ~31.9 GB download, WoLr2 ~94 GB, HUMAnN up to ~71 GB — see
> Section 4): they survive the cleanup and later runs reuse them. And
> because the work directory is deleted when a `low_disk` run succeeds, the run
> cannot be resumed afterwards: re-running with a changed or added option
> repeats the multi-hour eggNOG-mapper and read-level steps from scratch. The
> pipeline warns about this combination at launch; see
> [`DISK_OPTIMIZATION.md`](DISK_OPTIMIZATION.md).

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
│   ├── execution_trace_*.txt                   # Task trace: timing, %CPU, peak RSS, I/O, requested cpus/memory/time, attempt, work dir
│   ├── pipeline_dag_*.svg                      # Pipeline DAG
│   ├── software_versions.yml                   # Versions of every tool used in the run
│   └── contig_filtering_summary.txt            # Contig filtering summary (if assembly_mode != 'none')
├── 00_summary/                                 # Run-level summary (if results_summary=true)
│   ├── results_summary.html                    # Self-contained report: figures, sortable tables, methods
│   ├── results_summary.xlsx                    # Every table, plus a README sheet defining the columns
│   ├── README.md, README.docx                  # What the summary holds: stages run, tables, figures, units, caveats
│   ├── LLM_GUIDE.md                            # The summary folder (tables, figures, methods), for a language model or a script
│   ├── tables/                                 # samples_qc, taxonomy_long, taxonomy_<rank>_wide, taxonomy_species_prevalence,
│   │                                           # watchlist, arg_long, arg_gene_rpkm, arg_family_rpkm,
│   │                                           # arg_presence, arg_drug_class_rpkm, arg_priority_panel, group_arg_detection,
│   │                                           # arg_method_comparison, arg_method_overlap (K-MARVEL and RGI both ran),
│   │                                           # arg_sarg_long, arg_sarg_type_cpc, arg_sarg_subtype_cpc (ARGs-OAP),
│   │                                           # function_read_long, function_read_<ontology> (read functional)
│   └── figures/                                # qc_read_fate, taxonomy heatmaps and genus bars, watchlist, ARG heatmaps,
│                                               # arg_vs_depth, arg_method_upset, arg_method_by_category,
│                                               # arg_sarg_type_heatmap, function_cog_heatmap (.png/.pdf/.svg)
├── clean_reads/                                # Decontaminated reads (only if --store_clean_reads)
│   └── {sample}/                               # Per-sample clean R1/R2/Singleton FASTQs
├── 01_quality_control/                         # Quality control
│   ├── fastp/                                  # FastP reports per sample (if quality_control=true)
│   │   └── {sample}/                           # Per-sample QC results
│   │       ├── {sample}.fastp.html             # HTML report
│   │       ├── {sample}.fastp.json             # JSON report
│   │       └── {sample}.fastp.log              # Log file
│   └── summary/                                # Aggregated read counts (also with quality_control=false)
│       ├── Reads_report.csv                    # Read count summary
│       └── *.png                               # QC plots
├── 02_taxonomy/                                # Taxonomic profiling (if taxonomic_profiler != 'none')
│   ├── kraken2/                                # Kraken2 results (if taxonomic_profiler='kraken2')
│   │   └── {sample}_*.report.txt               # Per-sample Kraken2 reports
│   ├── bracken/                                # Bracken abundance estimates (kraken2 only)
│   │   ├── {sample}_{sample}.tsv               # Per-sample Bracken abundances
│   │   └── {sample}_{sample}.kraken2.report_bracken.txt  # Kraken-style report
│   ├── sourmash/                               # Sourmash results (if taxonomic_profiler='sourmash')
│   │   └── *.with-lineages.csv                 # Per-sample Sourmash results
│   ├── sylph/                                  # Sylph results (if taxonomic_profiler='sylph')
│   │   ├── {sample}_sylph_{db}.tsv             # Per-sample sylph genome profile (sylph profile -u)
│   │   └── {sample}_sylph_{db}.sylphmpa        # Per-sample lineage profile (sylph-tax, every rank down to t__ genome)
│   ├── tables/                                 # Taxonomy tables
│   │   ├── *.tsv                               # Abundance tables (OTU = relative abundance for sylph)
│   │   ├── sylph_{db}_read_counts.tsv          # Sylph: estimated read pairs per genome
│   │   └── Reads_report_taxonomy.csv           # Reads_report.csv plus the profiler's DB, Unclassified, Classified
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
│   │       ├── {sample}_contigs.fa             # Raw MEGAHIT contigs
│   │       └── {sample}_read_reduction.tsv     # Pairs in/kept, + {sample}_bbnorm.log (if assembly_read_reduction='normalize')
│   └── coassembly/                             # Co-assembly (if assembly_mode='coassembly')
│       ├── coassembly_filtered_contigs.fa      # Filtered contigs
│       ├── coassembly_contig.stats             # Assembly statistics
│       ├── coassembly_contigs.fa               # Raw MEGAHIT contigs
│       └── coassembly_read_reduction.tsv       # Pairs in/kept, + coassembly_bbnorm.log (if assembly_read_reduction='normalize')
├── 04_binning/                                 # Metagenomic binning (if include_binning=true)
│   ├── per_sample/                             # Per-sample binning (if assembly_mode='assembly')
│   │   └── {sample}/
│   │       ├── raw_bins/                       # Raw bins, one folder per selected binner
│   │       │   ├── metabat2/                   # MetaBAT2 bins
│   │       │   ├── semibin/                    # SemiBin bins
│   │       │   └── comebin/                    # COMEBin bins
│   │       ├── refined_bins/                   # Binette refined bins (if >= 2 binners)
│   │       ├── quality/                        # CheckM2 quality reports
│   │       │   └── checkm2/
│   │       └── taxonomy/                       # GTDB-TK taxonomy
│   │           └── gtdbtk/
│   ├── coassembly/                             # Co-assembly binning (if assembly_mode='coassembly')
│   │   ├── raw_bins/                           # Raw bins, one folder per selected binner
│   │   ├── refined_bins/                       # Binette refined bins (if >= 2 binners)
│   │   ├── quality/                            # CheckM2 quality reports
│   │   │   └── checkm2/
│   │   ├── taxonomy/                           # GTDB-TK taxonomy
│   │   │   └── gtdbtk/
│   │   ├── coverage/                           # Bin coverage per sample (<sample>_bin_depth.tsv)
│   │   └── summary/                            # Bin summary statistics (quality, taxonomy, coverage)
│   ├── quality/                                # Aggregated quality reports (if assembly_mode='assembly')
│   │   └── summary/
│   │       ├── *.csv                           # Quality summary tables
│   │       └── *.png                           # Quality plots
│   └── taxonomy/                               # Aggregated taxonomy reports (if assembly_mode='assembly')
│       └── summary/
│           ├── *.csv                           # Taxonomy summary tables
│           └── *.png                           # Taxonomy plots
├── 05_arg_prediction/                          # ARG predictions
│   ├── read_level/                             # Read-level ARG
│   │   ├── args_oap/                           # ARGs-OAP per sample (if args_oap_prediction=true)
│   │   │   └── {sample}_args_oap/              # metadata.txt (nRead, n16S, nCell), normalized_cell / normalized_16S / ppm tables by type, subtype, gene
│   │   ├── args_oap_summary/                   # ARGs-OAP run-level tables (if args_oap_prediction=true)
│   │   │   ├── args_oap_long.tsv               # One row per sample and SARG type, subtype or gene (copies per cell, per 16S, ppm)
│   │   │   └── args_oap_cells.tsv              # nRead, n16S, nCell, SARG reads and status per sample
│   │   ├── rgi/                                # RGI per-sample results (if rgi_prediction=true)
│   │   │   └── {sample}/
│   │   │       ├── {sample}_rgi_bwt.allele_mapping_data.txt
│   │   │       ├── {sample}_rgi_bwt.gene_mapping_data.txt
│   │   │       ├── {sample}_rgi_bwt.overall_mapping_stats.txt
│   │   │       ├── {sample}_rgi_bwt.reference_mapping_stats.txt
│   │   │       ├── {sample}_rgi_bwt.artifacts_mapping_stats.txt  # when RGI writes it
│   │   │       └── {sample}_rgi_bwt.sorted.length_100.bam        # absent for a sample without hits
│   │   ├── rgi_kmer/                           # RGI pathogen-of-origin (if rgi_prediction=true)
│   │   │   └── {sample}/
│   │   │       ├── *_61mer_analysis.json
│   │   │       ├── *_61mer_analysis.allele.txt
│   │   │       └── *_61mer_analysis.gene.txt
│   │   ├── rgi_summary/                        # RGI aggregated reports (if rgi_prediction=true)
│   │   │   ├── RGI_summary_report.csv
│   │   │   ├── RGI_amr_gene_family_distribution.png
│   │   │   ├── RGI_drug_class_profile.png
│   │   │   ├── RGI_resistance_mechanisms.png
│   │   │   └── RGI_sample_amr_counts.png
│   │   ├── kmarvel/                            # K-MARVEL per-sample results (if kmarvel_prediction=true)
│   │   │   └── {sample}/
│   │   │       ├── {sample}_kmarvel.tsv        # One row per assembled ARG (CARD ARO, drug class, identity, coverage, supporting read pairs, blast_* columns); header-only when none
│   │   │       ├── {sample}_kmarvel.json       # K-MARVEL JSON result
│   │   │       └── {sample}_kmarvel.log        # K-MARVEL stdout
│   │   └── kmarvel_summary/                    # K-MARVEL aggregated reports (if kmarvel_prediction=true)
│   │       ├── KMARVEL_summary_long.tsv        # One row per sample and ARG
│   │       ├── KMARVEL_supporting_reads_matrix.tsv  # ARG (CARD ARO) x sample supporting read pairs
│   │       ├── KMARVEL_gene_family_matrix.tsv  # AMR gene family x sample supporting read pairs (use this to compare samples)
│   │       └── KMARVEL_drug_class_matrix.tsv   # Drug class x sample number of ARGs
│   ├── contig_level/                           # Contig-level ARG (if contig_tax_and_arg=true)
│   │   ├── deeparg/                            # DeepARG predictions per sample
│   │   │   └── {sample}/                       # (coassembly/ under co-assembly)
│   │   │       └── *_contigs_deep_arg.out.mapping.ARG
│   │   ├── summary/                            # ARG summary reports
│   │   │   └── Contig_tax_and_arg_prediction.tsv
│   │   └── figures/                            # ARG visualization
│   │       └── *.png
│   └── bin_level/                              # Bin-level ARG (if arg_bin_clustering=true)
│       ├── proteins/                           # Per-bin Pyrodigal proteins (shared gene calling, split by bin)
│       │   └── {sample}/                       # (coassembly/ under co-assembly)
│       ├── deeparg/                            # Per-bin DeepARG hits (DEEPARG_CONTIGS split by bin)
│       │   └── {sample}/                       # (coassembly/ under co-assembly)
│       └── clustering/                         # MMseqs2 clustering results
│           └── *_cluster.tsv
├── 06_contig_taxonomy/                         # Contig taxonomy (if contig_tax_and_arg=true)
│   ├── metabuli/                               # Metabuli per contig
│   │   └── {sample}/                           # (coassembly/ under co-assembly)
│   │       ├── {sample}_classifications.tsv    # Per-contig call, score, rank and lineage
│   │       └── {sample}_report.tsv             # Clade report
│   ├── genomad/                                # geNomad (if genomad_prediction=true)
│   │   └── {sample}/                           # *_summary/ (plasmid/virus summaries, genes, sequences), *_aggregated_classification/, *_find_proviruses/, ...
│   └── figures/                                # BlobTools plots
│       └── *.png
└── 07_functional_annotation/                   # Functional annotation
    ├── gene_calling/                           # Pyrodigal ORF predictions on contigs
    │   └── {sample}/ or coassembly/            # (if contig_tax_and_arg, arg_bin_clustering or contig_level_functional)
    │       ├── {sample}.faa.gz
    │       ├── {sample}.fna.gz
    │       ├── {sample}.gff.gz
    │       └── {sample}.score.gz
    ├── eggnog/                                 # eggNOG-mapper v3 functional annotation
    │   └── {sample}/ or coassembly/            # (if contig_level_functional=true)
    │       ├── {sample}.emapper.seed_orthologs
    │       └── {sample}.emapper.annotations
    ├── dbcan/                                  # run_dbcan v5 CAZy annotation
    │   └── {sample}/ or coassembly/            # (if contig_level_functional=true and functional_cazy=true)
    │       ├── {sample}.overview.tsv           # per-gene CAZy calls with per-tool columns retained
    │       ├── {sample}.dbCANsub_hmm_results.tsv  # dbCAN-sub results incl. substrate predictions
    │       ├── {sample}.dbCAN_hmm_results.tsv  # raw dbCAN HMM results (provenance)
    │       └── {sample}.diamond.out            # raw DIAMOND-vs-CAZy results (provenance)
    ├── gene_abundance/                         # featureCounts per-gene read counts
    │   └── {sample}/                           # (if contig_level_functional=true; per sample in both assembly modes)
    │       ├── {sample}.featureCounts.txt      # Geneid, coordinates, Length, read count
    │       └── {sample}.featureCounts.txt.summary  # assigned vs unassigned alignments
    ├── microbecensus/                          # MicrobeCensus average genome size
    │   └── {sample}/                           # (if microbecensus=true and contig_level_functional=true or read_level_functional != none;
    │                                           #  per sample in both assembly modes)
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
    │   ├── read_function_wide_{ontology}_native.tsv  # wide native matrix (woltka: ko, ec, cog, pfam, metacyc read counts; superfocus: seed_level1/2/3 read counts; humann: ko, ec, metacyc RPK)
    │   ├── read_function_wide_{ontology}_cpge.tsv    # wide CPGE matrix (woltka and humann; blank columns for samples without AGS)
    │   ├── read_annotated_fraction.tsv         # share annotated, per ontology (woltka: of ORF-assigned reads; superfocus: of all reads; humann: reads mapped + RPK share per ontology)
    │   └── read_sample_summary.tsv             # reads assigned / ambiguous, AGS, CPGE status, tool + DB version
    ├── reads/                                  # read-level functional profiling, one backend per run
    │   ├── woltka/{sample}/                    # (if read_level_functional=woltka; needs no assembly)
    │   │   ├── {sample}.bowtie2.log            # alignment summary against WoLr2
    │   │   ├── {sample}.woltka_orf.tsv         # Woltka per-ORF read counts
    │   │   ├── {sample}.woltka_functions.tsv   # per-function read counts and RPK
    │   │   ├── {sample}.woltka_summary.tsv     # reads assigned / annotated per ontology
    │   │   └── {sample}.woltka_unassigned.tsv  # reads left unassigned by --woltka_uniq
    │   ├── superfocus/{sample}/                # (if read_level_functional=superfocus; needs no assembly)
    │   │   ├── {sample}.superfocus_functions.tsv   # SEED level 1-3 read counts (seed_level1/2/3)
    │   │   ├── {sample}.superfocus_summary.tsv     # reads given to SUPER-FOCUS / with an accepted hit
    │   │   ├── {sample}.superfocus_all_levels_and_function.xls  # raw SUPER-FOCUS table (tab-separated)
    │   │   ├── {sample}.superfocus_subsystem_level_{1,2,3}.xls  # raw SUPER-FOCUS per-level tables (not used, see note)
    │   │   └── {sample}.superfocus.log         # SUPER-FOCUS run log
    │   └── humann/{sample}/                    # (if read_level_functional=humann; needs no assembly)
    │       ├── {sample}_1_metaphlan_profile.tsv    # MetaPhlAn 4.1.2 prescreen profile
    │       ├── {sample}_2_genefamilies.tsv     # HUMAnN gene families (RPK; UniRef90 / UniClust90, stratified rows included)
    │       ├── {sample}_3_reactions.tsv        # HUMAnN MetaCyc reactions (RPK)
    │       ├── {sample}_4_pathabundance.tsv    # HUMAnN MetaCyc pathway abundance (RPK)
    │       ├── {sample}_0.log                  # HUMAnN run log
    │       ├── {sample}.humann_functions.tsv   # ko / ec / metacyc RPK composed by the pipeline
    │       └── {sample}.humann_summary.tsv     # reads given / mapped, RPK share annotated per ontology
    └── mags/                                   # Bakta MAG-level annotation, one file set per bin
        └── {sample}/ or coassembly/            # (if mag_level_functional=true; requires binning with >= 2 binners)
            ├── {sample}_{bin}.gff3             # annotation in GFF3 (coassembly_{bin}.* under co-assembly)
            ├── {sample}_{bin}.gbff             # annotation in GenBank flat file
            ├── {sample}_{bin}.faa              # protein sequences
            ├── {sample}_{bin}.fna              # replicon/contig sequences
            ├── {sample}_{bin}.tsv              # per-feature annotation table
            ├── {sample}_{bin}.txt              # per-bin annotation summary
            └── {sample}_{bin}.hypotheticals.tsv  # hypothetical-protein table (+ .hypotheticals.faa)
```

> **Gene identifiers and assembly mode:** in the gene abundance tables, `Geneid`
> is `<contig>_<n>`, matching the protein ids in `gene_calling/` and the query
> ids in the eggNOG annotations. Under `--assembly_mode coassembly` every
> sample is counted against the same co-assembly gene set, so gene ids are
> shared across samples and gene-level comparisons between samples within the
> run are valid. Under `--assembly_mode assembly` each sample has its own
> assembly and gene set: gene ids are sample-specific, and only function-level
> results (not per-gene rows) may be compared across samples. The directory
> keys follow the same split: under co-assembly, `gene_calling/`, `eggnog/`,
> `dbcan/` and `mags/` hold one set keyed `coassembly` (files named
> `coassembly.*` / `coassembly_<bin>.*`), while `gene_abundance/` and
> `microbecensus/` stay per sample in both modes. Counts are
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
> but over-predicts, while the assembly-based branch is more precise. The accession
> vocabularies do match the contig branch's, so the same term can be looked up in
> both: Woltka COG accessions are COG functional-category letters (the WoLr2 KO →
> COG ortholog ids are mapped with the same NCBI `cog-24.def.tab` table, which is
> downloaded for `--read_level_functional woltka` too; a COG with several categories
> contributes to each, and a read counts once per category even when its ORF
> reaches that category through several COGs). Nine WoLr2 `ko-to-cog` entries
> cannot be mapped (five malformed ids in the release, e.g. `COG:1140`, and four
> ids absent from COG 2024, e.g. `COG3632`); they are skipped with a warning in the
> task log, and the KO's other COGs still count. Pfam accessions are Pfam family
> names (`AAA`), with the versioned WoLr2 accession (`PF00004.32`) in
> `description`. MetaCyc accessions are pathways (`PWY-5101`); EC numbers are
> derived via KO.
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
>
> For HUMAnN, `abundance_native` is HUMAnN's RPK (`native_unit = rpk`;
> `--count-normalization RPKs`), unstratified: `metacyc` is the pathway
> abundance (accession = pathway id, description = its name; `UNMAPPED` /
> `UNINTEGRATED` excluded), `ko` / `ec` the summed RPK of the gene families
> carrying each term (`READS_UNMAPPED` excluded). `abundance_cpge` = RPK /
> genome equivalents, blank when MicrobeCensus is unavailable for the sample.

### Database Storage Directory

By default, databases are stored separately from results at `<output_dir>/../databases/`, one folder per database and registry key (the default keys are shown):

```
databases/                            # Database storage (configurable via --databases_dir)
├── bwamem3_index/phiX174_human/      # Combined bwa-mem3 decontamination index (host + phiX; _dense with --bwamem3_dense_sa)
├── bowtie_index/phiX174_human/       # Combined Bowtie2 index instead, under --read_aligner bowtie2
├── kraken/standard-8/                # Kraken2 database (--kraken2_db), k2_ref_db/
├── sourmash/gtdb_220_k31/            # sourmash taxonomy (taxonomy.sqldb)
├── sylph/gtdb_r232/                  # Sylph database (sylph_db/*.syl2db) + sylph-tax metadata (sylph_db/*.tsv.gz)
├── metabuli/refseq_standard/         # Metabuli database for contig taxonomy (metabuli_db/, ~110 GB, with its taxdump)
├── genomad/v1.9/                     # geNomad database 1.9 (genomad_db/, ~2.2 GB)
├── deeparg/1.0.4/                    # DeepARG database (deeparg_db/; unversioned upstream, keyed by the deeparg version)
├── rgi/latest_wildcard_k61/          # CARD database for RGI (card_db/; --rgi_card_version, _wildcard_k<--rgi_kmer_size> with WildCARD)
├── kmarvel/card_v4.0.1/              # CARD 4.0.1 protein models for K-MARVEL (kmarvel_db/, ~4.1 MB)
├── checkm2/v3/                       # CheckM2 database (*.dmnd)
├── gtdbtk/release_232/               # GTDB-TK database (gtdbtk_db/)
├── eggnog/emapper-3.0/               # eggNOG 7 data for eggNOG-mapper v3 (eggnog_db/, ~44 GB)
├── dbcan/db_v5-2-9_5-5-2026/         # dbCAN database for run_dbcan v5 (dbcan_db/, ~7.4 GB)
├── cog/cog-24/                       # NCBI COG definitions table (cog-24.def.tab, ~410 KB)
├── bakta/v6.0-full/                  # Bakta database v6.0 (bakta_db/; full 31.9 GB / light 1.3 GB download)
├── woltka/wolr2/                     # Web of Life WoLr2 subset for Woltka (woltka_db/, ~94 GB)
├── superfocus/db90/diamond/          # SUPER-FOCUS DB_90 root for the selected aligner (superfocus_db/, ~1.9-2.5 GB unpacked)
├── humann/v4_alpha-full/             # HUMAnN 4 database root + MetaPhlAn vOct22 (humann_db/, ~71 GB / ~33 GB EC-filtered)
└── stub/                             # Placeholder databases of -stub runs (same layout)
```

Copies of builds from your own files (not reused automatically, pass them with `--custom_*`): `bwamem3_index/contaminants_index` or `bowtie_index/contaminants_index` (from `--custom_host_fasta` / `--custom_phiX_fasta`), `rgi/card_db_with_wildcard` (from `--custom_rgi_wildcard`) and `sourmash/taxonomy.sqldb` (from a `--custom_sourmash_db` lineages file).

### Conditional Outputs

The following outputs are only generated when specific parameters are enabled:

| Output Directory | Required Parameter | Description |
|------------------|-------------------|-------------|
| `00_summary/` | `results_summary=true` (default) | Run-level summary: tables, Excel workbook, figures, HTML report |
| `01_quality_control/fastp/` | `quality_control=true` | FastP reports per sample |
| `01_quality_control/summary/` | always | Read-count report (`Reads_report.csv`; with a profiler, `02_taxonomy/tables/Reads_report_taxonomy.csv` adds the classified fraction) |
| `02_taxonomy/` | `taxonomic_profiler != 'none'` | Taxonomic profiling results |
| `02_taxonomy/phyloseq/*.RDS` | `create_phyloseq_rds=true` | R phyloseq object for downstream analysis |
| `03_assembly/` | `assembly_mode != 'none'` | Assembly results |
| `03_assembly/per_sample/` | `assembly_mode='assembly'` | Per-sample assemblies |
| `03_assembly/coassembly/` | `assembly_mode='coassembly'` | Co-assembly results |
| `04_binning/` | `include_binning=true` | Binning and bin refinement results |
| `05_arg_prediction/read_level/args_oap/`, `args_oap_summary/` | `args_oap_prediction=true` | ARGs-OAP (SARG) ARG profiles and read, 16S and cell estimates |
| `05_arg_prediction/read_level/rgi/` | `rgi_prediction=true` | RGI AMR predictions with pathogen-of-origin |
| `05_arg_prediction/read_level/kmarvel/` | `kmarvel_prediction=true` | K-MARVEL per-sample ARG tables and the `kmarvel_summary/` matrices |
| `05_arg_prediction/contig_level/` | `contig_tax_and_arg=true` | Contig-level ARG predictions |
| `05_arg_prediction/bin_level/` | `arg_bin_clustering=true` | Bin-level ARG clustering |
| `06_contig_taxonomy/` | `contig_tax_and_arg=true` | Contig taxonomic annotation (BlobTools) |
| `07_functional_annotation/gene_calling/` | `contig_tax_and_arg=true` or `contig_level_functional=true` | Pyrodigal ORF predictions on contigs (one `coassembly/` set under co-assembly) |
| `07_functional_annotation/eggnog/` | `contig_level_functional=true` | eggNOG-mapper functional annotation of predicted proteins (one `coassembly/` set under co-assembly) |
| `07_functional_annotation/dbcan/` | `contig_level_functional=true` and `functional_cazy=true` | run_dbcan CAZy annotation with per-tool calls and substrate predictions (one `coassembly/` set under co-assembly) |
| `07_functional_annotation/gene_abundance/` | `contig_level_functional=true` | featureCounts per-gene read counts (per sample in both assembly modes) |
| `07_functional_annotation/microbecensus/` | `microbecensus=true` and (`contig_level_functional=true` or `read_level_functional != 'none'`) | MicrobeCensus average genome size and genome equivalents (per sample in both assembly modes) |
| `07_functional_annotation/summary/` | `contig_level_functional=true` | Study-level tables: gene/function abundance with TPM and CPGE, wide matrices per ontology, annotated fraction, AGS summary |
| `07_functional_annotation/mags/` | `mag_level_functional=true` (needs `include_binning=true` and ≥2 binners) | Bakta per-bin MAG annotation (GFF3, GBFF, FAA, FNA, TSV, summary per refined bin; under `coassembly/` with co-assembly binning) |
| `07_functional_annotation/reads/woltka/` | `read_level_functional='woltka'` | Woltka read-level ORF and function tables per sample (no assembly needed) |
| `07_functional_annotation/reads/superfocus/` | `read_level_functional='superfocus'` | SUPER-FOCUS SEED level 1-3 tables and raw SUPER-FOCUS outputs per sample (no assembly needed) |
| `07_functional_annotation/reads/humann/` | `read_level_functional='humann'` | HUMAnN 4 raw tables (MetaPhlAn profile, gene families, reactions, pathway abundance) and the composed ko / ec / metacyc tables per sample (no assembly needed) |
| `07_functional_annotation/summary/read_*.tsv` | `read_level_functional != 'none'` | Study-level read-branch tables (same schema, `source=reads`), reported separately from the contig branch |

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

Complete analysis using the default single binner (no refinement):

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

### Example 3b: Full Pipeline with Binning + Binette Refinement

Run multiple binners and combine results with Binette for higher quality MAGs:

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
    --bin_min_completeness 50 \
    --bin_max_contamination 10 \
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
bin_min_completeness: 50
bin_max_contamination: 10
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
    --args_oap_prediction true \
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
    --custom_decontamination_index /shared/db/bwamem3_index/contaminants_index \
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
    --taxonomic_profiler kraken2 \
    --custom_kraken_db /shared/databases/bugbuster/kraken2_gtdb220 \
    --custom_checkm2_db /shared/databases/bugbuster/checkm2/uniref100.KO.1.dmnd \
    --custom_gtdbtk_db /shared/databases/bugbuster/release232 \
    --include_binning true \
    -profile singularity
```

### Example 13: Summarise several batches

Large studies are often run in batches (one results directory each). `bin/results_summary.py` reads the
published results only, so it merges finished batches without rerunning any tool. Sample names must be unique
across batches; taxa are matched by genome/taxon ID, and all batches must use the same taxonomy database.

The script needs Python 3 with pandas, numpy and matplotlib 3.7 or newer (its figure legends use a placement
older versions lack); openpyxl or xlsxwriter adds the Excel workbook, which is skipped with a warning without
them. The pipeline runs it with Python 3.13, pandas 3.1.0, numpy 2.5.4, matplotlib 3.11.2 and openpyxl 3.1.5
(`modules/local/results_summary/main.nf`; every report process uses this image), and it also runs with pandas 2. Run it from a BugBuster checkout: it imports the `summary_*.py`
modules next to it in `bin/` and reads its default tier rules and watchlist from `assets/`. Either way works:

**With the pipeline's container image**, which has everything installed:

```bash
IMG=community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:a1f6f67b03eef1d9
docker run --rm -u $(id -u):$(id -g) -v $PWD:$PWD -v /path/to/BugBuster:/bb:ro -v /dbs:/dbs:ro -w $PWD $IMG \
    python /bb/bin/results_summary.py \
        --results batch1/results batch2/results batch3/results \
        --outdir project_summary \
        --metadata samples.tsv --group-by subject --order-by visit \
        --kmarvel-db /dbs/kmarvel/card_v4.0.1 \
        --card-prevalence /dbs/rgi/wildcard/card_prevalence.txt
```

Under Singularity/Apptainer, use `singularity exec -B /path/to/BugBuster:/bb oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:cd89ac4e30745df7 python /bb/bin/results_summary.py ...`.

**With a local Python** that has the packages above (for example
`conda install -c conda-forge pandas numpy "matplotlib-base>=3.7" openpyxl`):

```bash
python /path/to/BugBuster/bin/results_summary.py \
    --results batch1/results batch2/results batch3/results \
    --outdir project_summary \
    --metadata samples.tsv --group-by subject --order-by visit \
    --kmarvel-db /dbs/kmarvel/card_v4.0.1 \
    --card-prevalence /dbs/rgi/wildcard/card_prevalence.txt
```

`--metadata`, `--group-by` and `--order-by` are optional. The `batch` column is each directory's name, or its
parent's when it is called `results` (`batch1` above); `--batch-labels` sets the labels instead, one per
`--results` directory in the same order. Without `--kmarvel-db`, mutation-dependent ARG hits
are reported as `not_checked`. RGI results (`05_arg_prediction/read_level/rgi/` and `rgi_kmer/`) are read when
present; `--rgi-min-coverage` sets their cutoff (default 80). ARGs-OAP (`05_arg_prediction/read_level/args_oap_summary/`),
the read functional profile (`07_functional_annotation/summary/read_*.tsv`) and MicrobeCensus
(`07_functional_annotation/microbecensus/`) are read when present; batches without them get `NA_not_run`. The
README and `LLM_GUIDE.md` say which stages ran in each batch. `results_summary.py --help` lists every option.

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
    custom_sylph_db          = '/shared/db/sylph/gtdb-r232-c200-dbv2.syl2db'
    custom_sylph_taxonomy    = '/shared/db/sylph/gtdb_r232_metadata.tsv.gz'
    custom_kraken_db         = '/shared/db/kraken2_standard8'   // used with --taxonomic_profiler kraken2
    custom_checkm2_db        = '/shared/db/checkm2/uniref100.KO.1.dmnd'
    custom_gtdbtk_db              = '/shared/db/release232'
    custom_decontamination_index  = '/shared/db/bwamem3_index/contaminants_index'
}

// The local executor's machine-wide pool is copied from max_cpus/max_memory
// before a -c file is read, so a -c file that raises them sets the pool too
executor {
    cpus   = 32
    memory = '256.GB'
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

> **Raising the resource limits.** `--max_cpus`, `--max_memory` and `--max_time` cap every
> task, and with the local executor `max_cpus`/`max_memory` are also the run's machine-wide
> pool. On the command line or in a `-params-file`, they set both. In a `-c` file, the
> `params` values cap each task, but the pool keeps the profile's values unless the file also
> has the `executor` block above. Without it, a task asking for more than the pool fails with
> `Process requirement exceeds available CPUs` (or `memory`). This matters most on top of
> `-profile test`, whose pool and caps are 2 CPUs / 6 GB / 6 h.

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
--taxonomic_profiler kraken2 --custom_kraken_db /path/to/manually/downloaded/db
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

*BugBuster v2.0.0dev - Built with Nextflow*
