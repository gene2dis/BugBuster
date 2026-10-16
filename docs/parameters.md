# BugBuster Parameter Reference

Complete reference for all BugBuster pipeline parameters.

---

## Table of Contents

1. [Input/Output Parameters](#inputoutput-parameters)
2. [Pipeline Execution Parameters](#pipeline-execution-parameters)
3. [Database Selection Parameters](#database-selection-parameters)
4. [Custom Database Paths](#custom-database-paths)
5. [Quality Control Parameters](#quality-control-parameters)
6. [Taxonomy Profiling Parameters](#taxonomy-profiling-parameters)
7. [Taxonomy Visualization Parameters](#taxonomy-visualization-parameters)
8. [Results Summary Parameters](#results-summary-parameters)
9. [Assembly Parameters](#assembly-parameters)
10. [Binning Parameters](#binning-parameters)
11. [ARG Prediction Parameters](#arg-prediction-parameters)
12. [Functional Annotation Parameters](#functional-annotation-parameters)
13. [Alignment Parameters](#alignment-parameters)
14. [Resource Limit Parameters](#resource-limit-parameters)
15. [Advanced Parameters](#advanced-parameters)

---

## Input/Output Parameters

### `--input`
- **Type**: String (file path)
- **Required**: Yes
- **Description**: Path to CSV samplesheet containing sample information
- **Format**: CSV file with columns: `sample`, `r1`, `r2`, `s` (optional)
- **Rules** (checked at launch): `sample` starts with a letter or digit and holds only letters, digits, `_`, `.` or `-` (it becomes file and directory names); sample ids are unique; read files end in `.fastq`, `.fq`, `.fastq.gz` or `.fq.gz`; `r1` and `r2` are different files, and `s` differs from both. Under `--read_level_functional humann` a sample id may not contain `s__` or `t__`. In co-assembly with `--quality_control false` the read file names must carry their role (see the [manual](manual.md#input-requirements))
- **Example**: `--input samplesheet.csv`

### `--output`
- **Type**: String (directory path)
- **Required**: Yes
- **Description**: Output directory where results will be saved
- **Example**: `--output ./results`

### `--publish_dir_mode`
- **Type**: String
- **Default**: `copy`
- **Options**: `copy`, `symlink`, `link`, `move`, `copyNoFollow`, `rellink`
- **Description**: Method used to save pipeline results to output directory
- **Example**: `--publish_dir_mode symlink`

---

## Pipeline Execution Parameters

### `--quality_control`
- **Type**: Boolean
- **Default**: `true`
- **Description**: Enable quality control and host read filtering steps
- **Example**: `--quality_control true`

### `--assembly_mode`
- **Type**: String
- **Default**: `assembly`
- **Options**: `assembly`, `coassembly`, `none`
- **Description**: Genome assembly strategy
  - `assembly`: Per-sample assembly
  - `coassembly`: Pool all samples for single assembly
  - `none`: Skip assembly
- **Example**: `--assembly_mode coassembly`

### `--taxonomic_profiler`
- **Type**: String
- **Default**: `sourmash`
- **Options**: `kraken2`, `sourmash`, `sylph`, `none`
- **Description**: Tool for taxonomic profiling at read level. `sylph` profiles the paired clean reads by containment ANI against a sylph genome database (`sylph profile -u`) and turns the genome hits into lineages with sylph-tax; like `kraken2`, it uses R1/R2 only (a singleton file would be a separate sylph sample)
- **Example**: `--taxonomic_profiler kraken2`

### `--include_binning`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable metagenomic binning and refinement steps
- **Example**: `--include_binning true`

### `--binners`
- **Type**: String (comma-separated)
- **Default**: `semibin`
- **Options**: Any combination of `comebin`, `semibin`, `metabat2` (at least one required)
- **Description**: Selects which binning tools to run. When ≥2 binners are specified, MetaWRAP bin refinement is automatically enabled. When only one binner is selected, MetaWRAP is skipped and the binner's output is used directly.
- **Examples**:
  - `--binners semibin` — single fast binner, no MetaWRAP (default)
  - `--binners semibin,metabat2` — two binners + MetaWRAP refinement
  - `--binners semibin,metabat2,comebin` — all three binners + MetaWRAP refinement
- **Note**: Only takes effect when `--include_binning true` is set

### `--read_arg_prediction`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable ARG and ARGV gene prediction at read level using KARGA and KARGVA
- **Example**: `--read_arg_prediction true`

### `--rgi_prediction`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable AMR gene prediction with pathogen-of-origin analysis using RGI and CARD database
- **Example**: `--rgi_prediction true`
- **Note**: See [`docs/RGI_WILDCARD_USAGE.md`](RGI_WILDCARD_USAGE.md) for manual database preparation

### `--kmarvel_prediction`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable read-level ARG detection with [K-MARVEL](https://bitbucket.org/amr-avengers/k-marvel): CARD protein ARGs are assembled from each sample's host-removed reads (R1, R2 and singletons pooled), chimeric assemblies are removed, variant models are kept only with a CARD resistance SNP, and every ARG is confirmed with blastp. Needs no assembly. Results in `05_arg_prediction/read_level/kmarvel/` and `kmarvel_summary/`
- **Example**: `--kmarvel_prediction true`
- **Note**: Runs from the pipeline-built image `ghcr.io/gene2dis/bugbuster-kmarvel` (no bioconda package). Independent of `--read_arg_prediction` and `--rgi_prediction`

### `--contig_tax_and_arg`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable contig-level taxonomy and ARG prediction using BlobTools and DeepARG. DeepARG runs on the proteins of the shared Pyrodigal gene-calling step, published to `07_functional_annotation/gene_calling/` — the same gene set `--contig_level_functional` uses (one gene-calling pass when both are on). The former `05_arg_prediction/contig_level/prodigal/` output no longer exists. The megablast search memory-maps the ~434 GB NCBI nt database and runs fastest when much of it fits in RAM; see [Contig taxonomy against NCBI nt](manual.md#contig-taxonomy-against-ncbi-nt)
- **Example**: `--contig_tax_and_arg true`

### `--contig_level_functional`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable the contig-level functional annotation branch: shared Pyrodigal gene calling on contigs (published to `07_functional_annotation/gene_calling/`), eggNOG-mapper v3 annotation of the predicted proteins (published to `07_functional_annotation/eggnog/`), run_dbcan v5 CAZy annotation of the same proteins (on by default, `--functional_cazy`; published to `07_functional_annotation/dbcan/`), per-gene abundance quantification with featureCounts over the gene coordinates (published to `07_functional_annotation/gene_abundance/`, per sample in both assembly modes; multi-mapping policy via `--featurecounts_multimap`), MicrobeCensus average genome size estimation for CPGE normalization (on by default, `--microbecensus`), and study-level aggregation into TPM and copies-per-genome-equivalent tables per functional ontology (KO, COG, EC, Pfam, CAZy — plus the dbCAN CAZy calls as a separate backend) with an annotated-fraction report and AGS summary (published to `07_functional_annotation/summary/`). Downloads the eggNOG 7 database (~44 GB) and the dbCAN database (~7.4 GB) on first use. Requires an assembly (`--assembly_mode assembly` or `coassembly`). Runs on every container engine: singularity/apptainer use the official eggNOG-mapper 3.0.0-beta6 Apptainer image, docker/podman and the cloud profiles use the pipeline-built `ghcr.io/gene2dis/bugbuster-eggnog-mapper` image of the same version (see the container note in the manual). Memory: the eggNOG-mapper search peaked at ~80 GB on a real 5-sample co-assembly — see [Functional annotation branches](manual.md#functional-annotation-branches) in the manual for observed per-step requirements
- **Example**: `--contig_level_functional true`

### `--mag_level_functional`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable MAG-level functional annotation: Bakta annotates every refined bin, one task per bin, publishing GFF3, GBFF, FAA, FNA, TSV and a summary per bin to `07_functional_annotation/mags/<sample>/` (files named `<sample>_<bin>.*`). Requires `--include_binning` **and at least two `--binners`** — MetaWRAP refinement and its completeness/contamination quality filter only run with ≥2 binners, and only quality-filtered bins are annotated (both requirements are validated at launch). Downloads the Bakta database on first use (`--bakta_db`: full 31.9 GB or light 1.3 GB download). Bins are expected to be bacterial; the pipeline does not detect or exclude archaeal/eukaryotic/viral bins — check the GTDB-Tk bin taxonomy report before interpreting annotations of non-bacterial bins. Independent of `--contig_level_functional`, and runs on any container engine (docker included)
- **Example**: `--mag_level_functional true --include_binning --binners semibin,metabat2`

### `--read_level_functional`
- **Type**: String
- **Default**: `none`
- **Options**: `woltka`, `superfocus`, `humann`, `none`
- **Description**: Optional read-level functional profiling of the host-removed reads, one backend per run (validated at launch). `woltka`: Bowtie2 alignment against the Web of Life release 2 genomes (WoLr2) with the SHOGUN multi-hit settings, Woltka ORF classification (a read is assigned to an ORF when ≥ 80 % of it lies inside; mates count separately; a read whose reported alignments — up to 16 within the score threshold — hit k ORFs counts 1/k toward each unless `--woltka_uniq`), then per-function read counts and RPK for KO, EC (via KO), COG functional categories (via KO; WoLr2's COG ids mapped with the `--cog_db` table, as in the contig branch), Pfam families (by name, as in the contig branch; the versioned accession is kept in `description`) and MetaCyc pathways — each read counted once toward each distinct term of its ORF. Per-sample tables go to `07_functional_annotation/reads/woltka/<sample>/`; study-level `read_function_abundance.tsv` (same schema as the contig branch, `source = reads`, native unit = reads, plus CPGE from MicrobeCensus), `read_function_wide_<ontology>_{native,cpge}.tsv`, `read_annotated_fraction.tsv` and `read_sample_summary.tsv` go to `07_functional_annotation/summary/`. Reported **separately** from the assembly-based tables, never merged (read-level profiling recovers the unassembled fraction but over-predicts). Needs no assembly (works with `--assembly_mode none`) and runs on any container engine. Downloads the WoLr2 subset (~94 GB, `--woltka_db`) on first use unless `--custom_woltka_db` is given; the alignment needs **≥ 68 GB RAM** per task (the WoLr2 index requirement)
  `superfocus`: SUPER-FOCUS 1.8 against its DB_90 SEED subsystem cluster database with DIAMOND blastx (default) or MMseqs2 (`--superfocus_aligner`). Each sample's R1, R2 and singleton reads are concatenated into one query (mates count separately); a read's equal-best-e-value hits passing the SUPER-FOCUS defaults (≥ 60 % identity, ≥ 15 aa, e-value ≤ 1e-5; `SUPERFOCUS` `ext.args`) are counted, each read with a hit contributing exactly 1, divided 1/k across the k distinct SEED (subsystem, function) assignments of its best hits. The pipeline sums these to SEED subsystem levels 1–3 (ontologies `seed_level1`, `seed_level2`, `seed_level3`, path-qualified accessions `L1` / `L1 | L2` / `L1 | L2 | L3`); SEED is never mapped to KO or EC. Per-sample tables go to `07_functional_annotation/reads/superfocus/<sample>/`; the same `read_*` summary tables are written, with only the three `read_function_wide_seed_level{1,2,3}_native.tsv` matrices — **no CPGE** (a SEED hit has no gene length, so no RPK; `cpge_status = not_applicable`, AGS/GE still reported). Downloads the DB_90 archive for the selected aligner only (~0.74 GB DIAMOND / ~0.9 GB MMseqs2, `--superfocus_db`) unless `--custom_superfocus_db` is given; runs on any container engine
  `humann`: HUMAnN 4.0.0a2 (an **alpha** release, pinned exactly) with its MetaPhlAn 4.1.2 prescreen against the MetaPhlAn database `mpa_vOct22_CHOCOPhlAnSGB_202403` (the only one this HUMAnN version accepts), nucleotide search against the prescreened ChocoPhlAn pangenomes, then translated search of the remaining reads against UniRef90 (EC-filtered, the only protein database HUMAnN 4 offers). Each sample's R1, R2 and singleton reads are concatenated into one input (mates count separately). HUMAnN runs with `--count-normalization RPKs`; the pipeline reports MetaCyc pathway abundance (`metacyc`) and regroups the gene families to KO (`ko`, UniRef90 families only — UniClust90 families never map to a KO) and level-4 EC (`ec`) itself, a family carrying several terms contributing its full RPK to each. KO tables are sparse: HUMAnN 4.0.0a2 ships a legacy UniRef90→KO map, and on the nf-core test reads only 17–18 % of the gene-family abundance mapped to a KO (EC 39–61 %); the per-ontology share is reported in `read_annotated_fraction.tsv`. `abundance_native` is RPK (`native_unit = rpk`) and CPGE = RPK / genome equivalents. Per-sample HUMAnN tables (`_1_metaphlan_profile`, `_2_genefamilies`, `_3_reactions`, `_4_pathabundance`, log) go to `07_functional_annotation/reads/humann/<sample>/`, with the same `read_*` summary tables (six wide matrices: ko/ec/metacyc × native/cpge). Sample names must not contain `s__` or `t__` (a HUMAnN 4.0.0a2 defect; rejected at launch). Downloads ~71 GB (`--humann_db v4_alpha-full`) or ~33 GB (`v4_alpha-ec_filtered`) unless `--custom_humann_db` is given; runs on any container engine; the MetaPhlAn prescreen loads a ~20 GB Bowtie2 index (HUMAnN peaked at ~23.5 GB per sample on the test data; see [Functional annotation branches](manual.md#functional-annotation-branches))
- **Note**: Independent of `--contig_level_functional` and `--mag_level_functional` — all three can be enabled in the same launch. MicrobeCensus then runs once per sample and feeds the CPGE of both the contig and the read tables, and the read tables stay in their own `read_*` files.
- **Example**: `--read_level_functional woltka --custom_woltka_db /shared/databases/wol2`; `--read_level_functional superfocus --custom_superfocus_db /shared/databases/superfocus`; `--read_level_functional humann --custom_humann_db /shared/databases/humann`

### `--arg_bin_clustering`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable ARG gene prediction and clustering for horizontal gene transfer inference
- **Example**: `--arg_bin_clustering true`

### `--min_read_sample`
- **Type**: Integer
- **Default**: `0`
- **Minimum**: `0`
- **Description**: Minimum number of reads required per sample after QC filtering. Samples below this threshold are excluded.
- **Example**: `--min_read_sample 10000000`

---

## Database Selection Parameters

### `--phiX_index`
- **Type**: String
- **Default**: `phiX174`
- **Options**: `phiX174` (the only built-in entry; use `--custom_phiX_fasta` for another genome)
- **Description**: PhiX genome selection for automatic download (genome FASTA; the pipeline builds the decontamination index for `--read_aligner`)
- **Size**: 5.4 kB
- **Example**: `--phiX_index phiX174`

### `--host_db`
- **Type**: String
- **Default**: `human`
- **Options**: `human` (T2T-CHM13v2.0, the only built-in entry; other hosts are passed with `--custom_host_fasta`)
- **Description**: Host genome for read filtering (genome FASTA; the pipeline builds the decontamination index for `--read_aligner`)
- **Size**: 940 MB (T2T-CHM13v2.0)
- **Example**: `--host_db human`

### `--kraken2_db`
- **Type**: String
- **Default**: `standard-8`
- **Options**: `standard-8`, `gtdb_220`
- **Description**: Kraken2 database selection
- **Size**: 7.5 GB (standard-8), 497 GB (gtdb_220)
- **Example**: `--kraken2_db gtdb_220`

### `--sourmash_db`
- **Type**: String
- **Default**: `gtdb_220_k31`
- **Description**: Sourmash database selection
- **Size**: 17 GB
- **Example**: `--sourmash_db gtdb_220_k31`

### `--sylph_db`
- **Type**: String
- **Default**: `gtdb_r232`
- **Options**: `gtdb_r232` (two-stage `.syl2db`, 19 GB), `gtdb_r226` (`.syldb`, 18 GB), `gtdb_r220` (`.syldb`, 13.1 GB)
- **Description**: Sylph database selection: GTDB representative genomes (c200 builds), downloaded together with the sylph-tax taxonomy metadata of the same release (~5-13 MB) into `<databases_dir>/sylph/`
- **Example**: `--sylph_db gtdb_r226`

### `--karga_db`
- **Type**: String
- **Default**: `megares`
- **Description**: KARGA reference database selection (MEGARes)
- **Size**: 9.2 MB
- **Example**: `--karga_db megares`

### `--kargva_db`
- **Type**: String
- **Default**: `kargva`
- **Description**: KARGVA reference database selection
- **Size**: 1.5 MB
- **Example**: `--kargva_db kargva`

### `--kmarvel_db`
- **Type**: String
- **Default**: `card_v4.0.1`
- **Options**: `card_v4.0.1`
- **Description**: CARD release for K-MARVEL: the CARD 4.0.1 data tarball (md5-verified), from which the protein homolog, overexpression and variant model FASTAs are concatenated into `CARD_protein_all.fasta`, with `aro_index.tsv` (as `metadata.tsv`) and `snps.txt`. CARD is free for academic use; commercial use needs a CARD licence
- **Size**: 4.6 MB download, 4.1 MB stored
- **Example**: `--kmarvel_db card_v4.0.1`

### `--blast_db`
- **Type**: String
- **Default**: `nt`
- **Description**: BLAST database selection for contig taxonomy
- **Size**: 434 GB
- **Example**: `--blast_db nt`

### `--taxdump_files`
- **Type**: String
- **Default**: `ncbi`
- **Description**: Taxonomy dump selection for BlobTools
- **Size**: 448 MB
- **Example**: `--taxdump_files ncbi`

### `--checkm2_db`
- **Type**: String
- **Default**: `v3`
- **Description**: CheckM2 database version
- **Size**: 2.9 GB
- **Example**: `--checkm2_db v3`

### `--gtdbtk_db`
- **Type**: String
- **Default**: `release_232`
- **Description**: GTDB-TK database release version. The pinned GTDB-Tk 2.7.2 accepts only GTDB R232 data (older R220/R226 packages require GTDB-Tk ≤ 2.6.1 and no longer work)
- **Size**: ~61 GB download
- **Example**: `--gtdbtk_db release_232`

### `--eggnog_db`
- **Type**: String
- **Default**: `emapper-3.0`
- **Description**: eggNOG 7 data selection for eggNOG-mapper v3 (contig-level functional annotation). The whole emapper 3.0.x series reuses this data directory
- **Size**: 44 GB (uncompressed)
- **Example**: `--eggnog_db emapper-3.0`

### `--dbcan_db`
- **Type**: String
- **Default**: `db_v5-2-9_5-5-2026`
- **Description**: dbCAN database release selection for run_dbcan v5 CAZy annotation (`--contig_level_functional` branch with `--functional_cazy`, the default). Downloads the four protein-mode CAZyme files (`CAZy.dmnd`, `dbCAN.hmm`, `dbCAN-sub.hmm`, `fam-substrate-mapping.tsv`) from the pinned dbCAN S3 release
- **Size**: 7.4 GB (uncompressed)
- **Example**: `--dbcan_db db_v5-2-9_5-5-2026`

### `--bakta_db`
- **Type**: String
- **Default**: `v6.0-full`
- **Options**: `v6.0-full`, `v6.0-light`
- **Description**: Bakta database flavor for MAG-level functional annotation (`--mag_level_functional`). Both are the pinned Zenodo v6.0 release (schema 6, required by Bakta 1.12.x). The light database has reduced annotation sources and **changes annotation results**: selecting it is always an explicit user choice and is recorded in provenance (the `bakta_db` entry in `pipeline_info/software_versions.yml`); no profile (including `low_disk`) ever switches it automatically
- **Size**: full 31.9 GB download / light 1.3 GB download
- **Example**: `--bakta_db v6.0-light`

### `--cog_db`
- **Type**: String
- **Default**: `cog-24`
- **Options**: `cog-24`
- **Description**: NCBI COG definitions table for the contig functional branch (`--contig_level_functional`) and the Woltka read backend (`--read_level_functional woltka`, which maps the COG ids of WoLr2's KO → COG map the same way, so `cog` means COG functional categories in both branches). eggNOG-mapper v3 / eggNOG 7 reports most genes' `COG_category` as a COG ortholog id (e.g. `COG1629`) rather than a functional category; the pipeline maps each id to its COG functional-category letters with this table (a COG with several categories contributes to each). Downloaded from `https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.def.tab` and verified against the release's `checksums.md5`; it covers every COG id of the eggNOG 7 database
- **Size**: ~410 KB
- **Example**: `--cog_db cog-24`

### `--woltka_db`
- **Type**: String
- **Default**: `wolr2`
- **Options**: `wolr2`
- **Description**: Web of Life release for the Woltka read-level backend (`--read_level_functional woltka`). Downloads, from the official public host (`https://ftp.microbio.me/pub/wol2`), only the files the backend reads, in the FTP layout: the Bowtie2 index (`databases/bowtie2/WoLr2.*.bt2l`, 93.6 GB), ORF coordinates and lengths (`proteins/`), and the KEGG / MetaCyc / Pfam function maps (`function/`). The published md5 checksums are verified; the pipeline bundles none of this data
- **Size**: ~94 GB
- **Example**: `--woltka_db wolr2`

### `--superfocus_db`
- **Type**: String
- **Default**: `db90`
- **Options**: `db90`
- **Description**: SUPER-FOCUS database for the SUPER-FOCUS read-level backend (`--read_level_functional superfocus`). Downloads the prebuilt DB_90 archive for the selected `--superfocus_aligner` only from figshare (open.flinders.edu.au, CC0): DIAMOND format 3 `90_clusters.db.dmnd` (zip ~0.74 GB, ~1.9 GB unpacked) or MMseqs2 `mmseqs_90.zip` (~0.9 GB, ~2.5 GB unpacked), plus `database_PKs.txt` (SEED subsystem levels per subsystem) from the SUPER-FOCUS v1.8 tag; both md5-verified. Stored as a database root at `<databases_dir>/superfocus/superfocus_db`. The cluster level is fixed at DB_90
- **Size**: ~0.74 GB (DIAMOND) / ~0.9 GB (MMseqs2) download
- **Example**: `--superfocus_db db90`

### `--humann_db`
- **Type**: String
- **Default**: `v4_alpha-full`
- **Options**: `v4_alpha-full`, `v4_alpha-ec_filtered`
- **Description**: HUMAnN database set for the HUMAnN read-level backend (`--read_level_functional humann`). Both keys download the HUMAnN v4_alpha UniRef90 EC-filtered DIAMOND database (0.94 GB), the v4 utility mapping (2.8 GB; KO / EC maps and the MetaCyc pathway files) and the MetaPhlAn database `mpa_vOct22_CHOCOPhlAnSGB_202403` with its prebuilt Bowtie2 index (~22.5 GB, md5-verified); they differ in the ChocoPhlAn pangenome database: `v4_alpha-full` (44.8 GB) or `v4_alpha-ec_filtered` (6.9 GB, only EC-annotated genes — fewer gene families and KOs). The choice is recorded in provenance (`read_sample_summary.tsv`, `software_versions.yml`). Stored as a database root at `<databases_dir>/humann/humann_db`; HUMAnN publishes no checksums for its own archives, so those are checked for layout only
- **Size**: ~71 GB (`v4_alpha-full`) / ~33 GB (`v4_alpha-ec_filtered`) download
- **Example**: `--humann_db v4_alpha-ec_filtered`

### `--databases_dir`
- **Type**: String (directory path)
- **Default**: `<output>/../databases`
- **Description**: Directory for storing downloaded databases (separate from results output)
- **Example**: `--databases_dir /shared/databases/bugbuster`

---

## Output & Cleanup Options

### `--store_clean_reads`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Publish decontaminated reads to `<output>/clean_reads/<sample>/`. This is publishing only, not caching — `-resume` remains the mechanism for reusing completed work.
- **Example**: `--store_clean_reads true`

> **Note**: Automatic work-dir cleanup is not a pipeline parameter — it comes from Nextflow's `cleanup = true` setting, enabled by `-profile low_disk` (or a custom config). See [`DISK_OPTIMIZATION.md`](DISK_OPTIMIZATION.md).

---

## Custom Database Paths

Override automatic downloads by providing custom database paths:

### `--custom_decontamination_index`
- **Type**: String (directory path)
- **Description**: Path to a pre-built combined decontamination index directory (host + phiX) matching `--read_aligner`: a bwa-mem3/bwa-mem2 index (`*.bwt.2bit.64`) under the default `bwamem3`, or a Bowtie2 index under `--read_aligner bowtie2`. A Bowtie2 index under `bwamem3` is rejected at launch: add `--read_aligner bowtie2` to keep using an index built by an earlier release. Skips both the genome downloads and the index build. Pointing it at the index the pipeline published (`<databases_dir>/bwamem3_index/contaminants_index`, or `<databases_dir>/bowtie_index/contaminants_index` under `bowtie2`) avoids rebuilding it in later runs.
- **Example**: `--custom_decontamination_index /path/to/bwamem3_index/contaminants_index`
- **Note**: See [`docs/DECONTAMINATION_QUICK_REFERENCE.md`](DECONTAMINATION_QUICK_REFERENCE.md) for details

### `--custom_phiX_fasta`
- **Type**: String (file path)
- **Description**: Path to a custom phiX genome FASTA; the pipeline builds the decontamination index (`--read_aligner`) from it
- **Example**: `--custom_phiX_fasta /path/to/phiX174.fasta`

### `--custom_host_fasta`
- **Type**: String (file path)
- **Description**: Path to a custom host genome FASTA, plain or gzipped (e.g. a non-human host), or a comma-separated list of them (spaces around the commas are ignored; paths may not contain commas or spaces). The pipeline builds one decontamination index (`--read_aligner`) from them and the phiX genome. The bwa-mem3 build requests ~24 GB RAM per Gbp of genome, capped by `--max_memory`, and the launch memory warning covers only the built-in `--host_db`, so size `--max_memory` yourself for large or several genomes. See [`YAML_PARAMETERS_GUIDE.md`](YAML_PARAMETERS_GUIDE.md)
- **Example**: `--custom_host_fasta /path/to/mouse_genome.fa` or `--custom_host_fasta /path/to/human.fa.gz,/path/to/mouse.fa.gz`

> **Removed**: the old `--custom_phiX_index` and `--custom_bowtie_host_index` parameters no longer exist (passing them aborts the run) — use `--custom_decontamination_index`, `--custom_phiX_fasta`, or `--custom_host_fasta` instead.

### `--custom_kraken_db`
- **Type**: String (directory path)
- **Description**: Path to custom Kraken2 database directory
- **Example**: `--custom_kraken_db /path/to/kraken2_db`

### `--custom_sourmash_db`
- **Type**: Array of strings
- **Description**: Paths to custom Sourmash database files [kmer_file, lineages_file]
- **Example**: `--custom_sourmash_db '["/path/to/kmers.zip", "/path/to/lineages.csv"]'`

### `--custom_sylph_db`
- **Type**: String (file path)
- **Description**: Path to a local sylph genome database (`.syldb` or `.syl2db`, e.g. a [pre-built database](https://sylph-docs.github.io/pre%E2%80%90built-databases/) or one made with `sylph sketch`). Must be given together with `--custom_sylph_taxonomy` (checked at launch), so a custom database is never paired with the registry's taxonomy; the two skip both downloads. Outputs are then labelled with the database's file name without extension (e.g. `sylph_gtdb-r232-c200-dbv2_otu_table.tsv` and `gtdb-r232-c200-dbv2` in the `Sylph DB` column) instead of the `--sylph_db` key
- **Example**: `--custom_sylph_db /path/to/gtdb-r232-c200-dbv2.syl2db`

### `--custom_sylph_taxonomy`
- **Type**: String (file path)
- **Description**: Path to a local sylph-tax metadata file matching the sylph database: genome accession, a tab, and a `d__...;p__...;...;s__...` lineage per line (`.tsv` or `.tsv.gz`). Use the file `sylph-tax download` fetches for a pre-built database (e.g. `gtdb_r232_metadata.tsv.gz`), or write your own for custom genomes. Must be given together with `--custom_sylph_db`
- **Example**: `--custom_sylph_taxonomy /path/to/gtdb_r232_metadata.tsv.gz`

### `--custom_checkm2_db`
- **Type**: String (file path)
- **Description**: Path to custom CheckM2 database file
- **Example**: `--custom_checkm2_db /path/to/checkm2.dmnd`

### `--custom_gtdbtk_db`
- **Type**: String (directory path)
- **Description**: Path to the directory that **directly contains** the unarchived GTDB-Tk reference data (`markers/`, `skani/`, `taxonomy/`, `msa/`, ...) — the pipeline sets `GTDBTK_DATA_PATH` to this directory. For the official packages that is the extracted release directory itself (e.g. `release232/`), not its parent. The pipeline pins GTDB-Tk 2.7.2, which per the upstream compatibility table accepts **only GTDB R232** data; older R220/R226 packages require GTDB-Tk ≤ 2.6.1 and will not work
- **Example**: `--custom_gtdbtk_db /path/to/release232`

### `--custom_deeparg_db`
- **Type**: String (directory path)
- **Description**: Path to custom DeepARG database directory
- **Example**: `--custom_deeparg_db /path/to/deeparg_db`

### `--custom_blast_db`
- **Type**: String (directory path)
- **Description**: Path to custom BLAST NT database directory
- **Example**: `--custom_blast_db /path/to/blast_nt`

### `--custom_taxdump_files`
- **Type**: String (directory path)
- **Description**: Path to custom NCBI taxdump directory
- **Example**: `--custom_taxdump_files /path/to/taxdump`

### `--custom_karga_db`
- **Type**: String (file path)
- **Description**: Path to custom KARGA database FASTA file
- **Example**: `--custom_karga_db /path/to/megares.fasta`

### `--custom_kargva_db`
- **Type**: String (file path)
- **Description**: Path to custom KARGVA database FASTA file
- **Example**: `--custom_kargva_db /path/to/kargva.fasta`

### `--custom_rgi_card_db`
- **Type**: String (directory path)
- **Description**: Path to pre-prepared CARD database directory for RGI: its `localDB/` (created by `rgi load --local`) must include WildCARD and the k-mer database for `--rgi_kmer_size` (`localDB/61mer_database.json`, `localDB/amr_61mer.txt`), which the pathogen-of-origin step needs; that step fails when they are missing. KMA indices are built by RGI when absent. Used alone, the directory is linked, not copied
- **Example**: `--custom_rgi_card_db /path/to/card_database`
- **Note**: See [`docs/RGI_WILDCARD_USAGE.md`](RGI_WILDCARD_USAGE.md) for database preparation instructions

### `--custom_rgi_wildcard`
- **Type**: String (directory path)
- **Description**: Path to WildCARD directory containing variant files. Use together with `--custom_rgi_card_db` to combine separate databases
- **Example**: `--custom_rgi_wildcard /path/to/wildcard_directory`
- **Requirements**: Directory must contain `index-for-model-sequences.txt`, the variant FASTA files and the k-mer database for `--rgi_kmer_size` (`61_kmer_db.json`, `all_amr_61mers.txt`), all decompressed. Both directories are copied into the preparation task and the combined database is published as a full copy, so loading the database once by hand and passing `--custom_rgi_card_db` alone saves that space
- **Note**: See [`docs/RGI_WILDCARD_USAGE.md`](RGI_WILDCARD_USAGE.md) for detailed usage examples

### `--custom_kmarvel_db`
- **Type**: String (directory path)
- **Description**: Path to a K-MARVEL database directory; skips the download. A database published by an earlier run (`databases/kmarvel/kmarvel_db`) can be reused as is
- **Example**: `--custom_kmarvel_db /path/to/databases/kmarvel/kmarvel_db`
- **Requirements**: Directory must contain `CARD_protein_all.fasta` (CARD protein homolog + overexpression + variant model FASTAs), `metadata.tsv` (CARD `aro_index.tsv`) and `snps.txt` (CARD `snps.txt`). No BLAST index is needed: each task indexes its own copy of the FASTA. An optional `DB_VERSION` file (its first line) is recorded as `card:` under `KMARVEL` in `pipeline_info/software_versions.yml`; without it the entry is `custom (CARD_protein_all.fasta md5 <md5>)`, so the reference can still be traced

### `--custom_dbcan_db`
- **Type**: String (directory path)
- **Description**: Path to custom dbCAN database directory in the run_dbcan v5 layout: `CAZy.dmnd`, `dbCAN.hmm`, `dbCAN-sub.hmm` (hyphen — the tool's expected filename), `fam-substrate-mapping.tsv`. An optional `DB_VERSION` file (one line, the release string) feeds provenance; without it the recorded database version is `custom`
- **Example**: `--custom_dbcan_db /path/to/dbcan_db`

### `--custom_eggnog_db`
- **Type**: String (directory path)
- **Description**: Path to custom eggNOG 7 data directory in the emapper-3.0 layout: `eggnog.db` (plus its `.fieldpresence.bin` and `.taxids.bin` caches), `eggnog.taxa.db` (plus `.traverse.pkl`), `eggnog_proteins.dmnd`, `go-basic.obo`. Must be eggNOG 7 data — eggNOG-mapper v3 rejects eggNOG 5 databases
- **Example**: `--custom_eggnog_db /path/to/emapper-3.0/data`

### `--custom_bakta_db`
- **Type**: String (directory path)
- **Description**: Path to custom Bakta database directory in the schema 6 layout (the content of an extracted `db.tar.xz`/`db-light.tar.xz`: `version.json`, `amrfinderplus-db/`, and the Bakta annotation databases). Must be schema 6 — Bakta 1.12.x rejects older schemas. An optional `DB_VERSION` file (one line) feeds provenance; without it the recorded database version is `custom`. **Note:** the official v6.0 tarball bundles an AMRFinderPlus database too old for the AMRFinderPlus in the pinned bakta container — refresh it once with `amrfinder_update --force_update --database <db>/amrfinderplus-db` (run inside the bakta container; see `docs/troubleshooting.md`), otherwise every annotation fails at the AMR expert step. The pipeline's auto-download path does this refresh automatically
- **Example**: `--custom_bakta_db /path/to/bakta_db`

### `--custom_cog_db`
- **Type**: String (file path)
- **Description**: Path to a local NCBI COG definitions table in the `cog-24.def.tab` layout (tab-separated, no header: COG id, functional-category letters, name, ...). It must cover every COG id the eggNOG database reports — an id missing from the table stops the aggregation with an error naming it — and, with `--read_level_functional woltka`, every COG id of the WoLr2 `function/kegg/ko-to-cog.map` (apart from nine known WoLr2 defects that are skipped; any other unmapped id stops `WOLTKA_CLASSIFY`)
- **Example**: `--custom_cog_db /shared/databases/cog/cog-24.def.tab`

### `--custom_woltka_db`
- **Type**: String (directory path)
- **Description**: Path to a local Web of Life release 2 mirror in the FTP layout of `https://ftp.microbio.me/pub/wol2` — at least `databases/bowtie2/WoLr2.{1,2,3,4,rev.1,rev.2}.bt2l`, `proteins/coords.txt.xz`, `proteins/length.map.xz`, `function/kegg/{orf-to-ko.map.xz,ko-to-ec.map,ko-to-cog.map,ko_name.txt}`, `function/metacyc/{orf-to-protein.map.xz,protein-to-enzrxn.map,enzrxn-to-reaction.map,reaction-to-pathway.map,pathway_name.txt}` and `function/pfam/{orf-to-pfam.map.xz,pfam_name.txt}` (a download recipe is in `docs/manual.md`, Manual Database Download). The index and the coordinate/map files must come from the same release. An optional `DB_VERSION` file (one line) feeds provenance; without it the recorded version is `custom (<index name>)`
- **Example**: `--custom_woltka_db /shared/databases/wol2`

### `--custom_superfocus_db`
- **Type**: String (directory path)
- **Description**: Path to a local SUPER-FOCUS database **root** — the directory that CONTAINS `db/`: `db/database_PKs.txt` plus `db/static/diamond/90_clusters.db.dmnd` (for `--superfocus_aligner diamond`) or `db/static/mmseqs2/90_clusters.db*` (for `mmseqs2`). Pointing at the `db/` folder itself fails with an explicit message, as does a database without the selected aligner's DB_90 files. A database root from an existing SUPER-FOCUS install (e.g. made with an older SUPER-FOCUS version; DIAMOND format 3 `.dmnd`) works unmodified. A download recipe is in `docs/manual.md`, Manual Database Download. An optional `DB_VERSION` file (one line) feeds provenance; without it the recorded version is `custom (<aligner> DB_90)`
- **Example**: `--custom_superfocus_db /shared/databases/superfocus`

### `--custom_humann_db`
- **Type**: String (directory path)
- **Description**: Path to a local HUMAnN database **root** containing `chocophlan/` (ChocoPhlAn v4_alpha pangenomes), `uniref/` (the UniRef90 EC-filtered `.dmnd`), `utility_mapping/` (the `full_mapping_v4_alpha` files: `map_ko_uniref90.txt.gz`, `map_level4ec_uniclust90.txt.gz`, their name maps and the two MetaCyc pathway files) and `metaphlan/` (`mpa_vOct22_CHOCOPhlAnSGB_202403.pkl` plus its `.bt2l` Bowtie2 index — HUMAnN 4.0.0a2 refuses any other MetaPhlAn database, e.g. vJun23). This is the layout `humann_databases --download` produces for the three HUMAnN archives, plus the MetaPhlAn files in their own folder. An optional `DB_VERSION` file (one line) feeds provenance; without it the recorded version is `custom`
- **Example**: `--custom_humann_db /shared/databases/humann`

---

## Quality Control Parameters

### `--fastp_n_base_limit`
- **Type**: Integer
- **Default**: `5`
- **Description**: Maximum number of N bases allowed in a read
- **Example**: `--fastp_n_base_limit 10`

### `--fastp_unqualified_percent_limit`
- **Type**: Integer
- **Default**: `10`
- **Range**: 0-100
- **Description**: Maximum percentage of unqualified bases allowed
- **Example**: `--fastp_unqualified_percent_limit 15`

### `--fastp_qualified_quality_phred`
- **Type**: Integer
- **Default**: `20`
- **Description**: Phred quality score threshold for qualified bases
- **Example**: `--fastp_qualified_quality_phred 25`

### `--fastp_cut_front_window_size`
- **Type**: Integer
- **Default**: `4`
- **Description**: Window size for cutting from front of reads
- **Example**: `--fastp_cut_front_window_size 5`

### `--fastp_cut_front_mean_quality`
- **Type**: Integer
- **Default**: `20`
- **Description**: Mean quality threshold for front cutting
- **Example**: `--fastp_cut_front_mean_quality 25`

### `--fastp_cut_right_window_size`
- **Type**: Integer
- **Default**: `4`
- **Description**: Window size for cutting from right of reads
- **Example**: `--fastp_cut_right_window_size 5`

### `--fastp_cut_right_mean_quality`
- **Type**: Integer
- **Default**: `20`
- **Description**: Mean quality threshold for right cutting
- **Example**: `--fastp_cut_right_mean_quality 25`

---

## Taxonomy Profiling Parameters

### `--kraken_confidence`
- **Type**: Number
- **Default**: `0.1`
- **Range**: 0-1
- **Description**: Kraken2 confidence threshold for taxonomic assignment
- **Example**: `--kraken_confidence 0.2`

### `--bracken_read_len`
- **Type**: Integer
- **Default**: `150`
- **Description**: Read length for Bracken abundance estimation
- **Example**: `--bracken_read_len 100`

### `--bracken_tax_level`
- **Type**: String
- **Default**: `S`
- **Options**: `D`, `P`, `C`, `O`, `F`, `G`, `S`
- **Description**: Taxonomic level for Bracken (D=Domain, P=Phylum, C=Class, O=Order, F=Family, G=Genus, S=Species)
- **Example**: `--bracken_tax_level G`

### `--sourmash_tax_rank`
- **Type**: String
- **Default**: `species`
- **Options**: `genus`, `species`, `strain`
- **Description**: Taxonomic rank for Sourmash classification
- **Example**: `--sourmash_tax_rank genus`

### `--sylph_read_seq_id`
- **Type**: Number (percent, > 0 and <= 100)
- **Default**: unset (sylph's automatic estimate)
- **Description**: Read sequence identity passed to `sylph profile --read-seq-id`. It drives sylph's unknown-fraction estimate (`-u`), i.e. the `Classified` column of `Reads_report.csv` and the estimated read pairs per genome; relative abundance is not affected. Reads from genomes absent from the database can skew the automatic estimate towards 100% classified (on simulated reads with 1/3 unknown: automatic 100%, `100` gave 66.1% against a true 66.7%), so set it when `Classified` looks too high
- **Example**: `--sylph_read_seq_id 99.5` (typical Illumina)

### Sylph outputs
Sylph runs with sylph's defaults plus `-u` (estimate the unknown fraction) and
`--read-seq-id` when `--sylph_read_seq_id` is set. Per sample it publishes the genome
profile and the sylph-tax lineage profile under `02_taxonomy/sylph/`. The
`02_taxonomy/tables/` set holds:

- `sylph_<db>_otu_table.tsv`: sylph's **relative (taxonomic) abundance** of each
  genome (`t__` accession), coverage-normalised, MetaPhlAn-like; this is the
  table behind the plots and the phyloseq objects
- `sylph_<db>_read_counts.tsv`: **estimated read pairs per genome**,
  `sequence_abundance / 100 × the read pairs given to sylph`. With `-u`,
  `sequence_abundance` is the genome's share of all input reads (unknown
  included), so the column sums to the classified pairs
- The `Classified` column of `Reads_report.csv` is the sum of `sequence_abundance / 100`
  (capped at 1.0); see `--sylph_read_seq_id` when it looks too high

---

## Taxonomy Visualization Parameters

### `--taxonomy_plot_levels`
- **Type**: String
- **Default**: `Phylum,Family,Genus,Species`
- **Description**: Comma-separated list of taxonomic levels to plot
- **Example**: `--taxonomy_plot_levels "Phylum,Class,Order,Family"`

### `--taxonomy_top_n_taxa`
- **Type**: Integer
- **Default**: `10`
- **Minimum**: `1`
- **Description**: Maximum taxa per plot. Each sample's most abundant taxa are included first, then the most abundant overall, so a taxon that dominates a few samples is not hidden in "Other"
- **Example**: `--taxonomy_top_n_taxa 20`

### `--create_phyloseq_rds`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Generate R phyloseq RDS files for downstream analysis
- **Example**: `--create_phyloseq_rds true`

---

## Results Summary Parameters

The run-level summary (`00_summary/`) consolidates the QC, taxonomy and K-MARVEL results into TSV tables, an
Excel workbook, figures and a self-contained HTML report. Counts are read pairs, ARG abundance is RPM/RPKM of
clean read pairs, and ARG calls carry a reporting tier. See the manual, section 6.16, for how the results are reported.

### `--results_summary`
- **Type**: Boolean
- **Default**: `true`
- **Description**: Write the run-level summary to `00_summary/`
- **Example**: `--results_summary false`

### `--sample_metadata`
- **Type**: String (path to a `.tsv`, `.csv` or `.txt` file)
- **Default**: none
- **Description**: Table with a `sample` column (or sample names in the first column) and any metadata columns; the columns are carried into every summary table
- **Example**: `--sample_metadata samples.tsv`

### `--summary_group_by`
- **Type**: String
- **Default**: none
- **Description**: Metadata column that groups samples (e.g. subject, site) in the summary figures and in `group_arg_detection.tsv`. Needs `--sample_metadata`
- **Example**: `--summary_group_by subject`

### `--summary_order_by`
- **Type**: String
- **Default**: none
- **Description**: Metadata column that orders samples within a group (e.g. visit, day; `V2` sorts before `V10`). With `--summary_group_by` it orders the samples of each group and the per-group tier 1 detection patterns (`group_arg_detection`). Needs `--sample_metadata`
- **Example**: `--summary_order_by visit`

### `--summary_batch_label`
- **Type**: String (letters, digits, `.`, `_`, `-`)
- **Default**: the name of the `--output` directory
- **Description**: Label written in the `batch` column of the summary tables. Set it when the batches of a project share an output directory name (e.g. every batch writes to `results/`), so that merged summaries tell the batches apart
- **Example**: `--summary_batch_label batch1`

### `--summary_min_clean_pairs`
- **Type**: Integer
- **Default**: `1000000`
- **Minimum**: `0`
- **Description**: Samples with fewer clean read pairs are flagged as low depth. A gene or taxon not found in such a sample is reported as `NA_low_depth`, not as absent. Samples are never dropped
- **Example**: `--summary_min_clean_pairs 5000000`

### `--summary_top_n`
- **Type**: Integer
- **Default**: `40`
- **Minimum**: `1`
- **Description**: Maximum rows in the summary ARG heatmap (the taxonomy heatmaps use `--summary_taxa_per_sample`)
- **Example**: `--summary_top_n 60`

### `--summary_taxa_per_sample`
- **Type**: Integer
- **Default**: `5`
- **Minimum**: `1`
- **Description**: The summary species and genus heatmaps show the union of each sample's N most abundant taxa, with no row cap, grouped by phylum and alphabetical within each phylum
- **Example**: `--summary_taxa_per_sample 10`

### `--summary_depth_marks`
- **Type**: String (comma-separated read-pair counts)
- **Default**: `10000000,20000000`
- **Description**: Clean read-pair depths drawn as dotted reference lines in the summary QC figure, counted in the report header ("samples with >=10M clean pairs") and used for the `depth_class` column. `''` for none. The low-depth threshold (`--summary_min_clean_pairs`) is always drawn
- **Example**: `--summary_depth_marks 5000000,10000000,20000000`

### `--card_prevalence`
- **Type**: String (path)
- **Default**: none
- **Description**: WildCARD `card_prevalence.txt` (in an RGI WildCARD database: `wildcard/card_prevalence.txt`). Adds to each ARG call the species of the same sample that CARD finds carrying that ARO, with their abundance, and whether CARD finds the ARO on plasmids. This is co-occurrence, not linkage
- **Example**: `--card_prevalence /dbs/rgi/wildcard/card_prevalence.txt`

### `--arg_priority_table`
- **Type**: String (path)
- **Default**: `assets/arg_priority_tiers.tsv`
- **Description**: ARG reporting tier rules: regular expressions on CARD ARO names or gene families, first match wins. Copy and edit the default to change the tiers
- **Example**: `--arg_priority_table my_tiers.tsv`

### `--summary_watchlist`
- **Type**: String (path)
- **Default**: `assets/summary_watchlist.tsv`
- **Description**: Taxa (species or genus) reported in the summary watchlist table and figure; GTDB suffixes are ignored when matching
- **Example**: `--summary_watchlist my_taxa.tsv`

---

## Assembly Parameters

### `--bbmap_length`
- **Type**: Integer
- **Default**: `1000`
- **Minimum**: `0`
- **Description**: Minimum contig length after BBMap filtering
- **Example**: `--bbmap_length 1500`

---

## Binning Parameters

### Basic Binning Parameters

#### `--binners`
- **Type**: String (comma-separated list)
- **Default**: `semibin`
- **Options**: Any combination of `comebin`, `semibin`, `metabat2` (at least one required)
- **Description**: Controls which binning tools are run. When ≥2 binners are selected, MetaWRAP bin refinement is automatically enabled. When only one binner is selected, its output is used directly (no MetaWRAP). Only takes effect when `--include_binning true`.
- **CLI examples**:
  ```bash
  --binners semibin                    # single binner, no MetaWRAP (default)
  --binners semibin,metabat2           # two binners + MetaWRAP refinement
  --binners semibin,metabat2,comebin   # all three binners + MetaWRAP refinement
  ```
- **YAML example**:
  ```yaml
  include_binning: true
  binners: "semibin,metabat2"
  ```

#### `--metabat_minContig`
- **Type**: Integer
- **Default**: `2500`
- **Description**: Minimum contig length for MetaBAT2 binning
- **Example**: `--metabat_minContig 3000`

#### `--metawrap_completeness`
- **Type**: Integer
- **Default**: `50`
- **Range**: 0-100
- **Description**: Minimum bin completeness threshold for MetaWRAP refinement
- **Example**: `--metawrap_completeness 70`

#### `--metawrap_contamination`
- **Type**: Integer
- **Default**: `10`
- **Range**: 0-100
- **Description**: Maximum bin contamination threshold for MetaWRAP refinement
- **Example**: `--metawrap_contamination 5`

#### `--semibin_env_model`
- **Type**: String
- **Default**: `human_gut`
- **Options**: `human_gut`, `dog_gut`, `ocean`, `soil`, `cat_gut`, `human_oral`, `mouse_gut`, `pig_gut`, `built_environment`, `wastewater`, `chicken_caecum`, `global`
- **Description**: SemiBin environment model for binning
- **Example**: `--semibin_env_model ocean`

### Advanced MetaBAT2 Parameters

#### `--metabat_maxP`
- **Type**: Integer
- **Default**: `95`
- **Range**: 0-100
- **Description**: Maximum percentage of good contigs for MetaBAT2
- **Example**: `--metabat_maxP 90`

#### `--metabat_minS`
- **Type**: Integer
- **Default**: `60`
- **Description**: Minimum score for MetaBAT2 binning
- **Example**: `--metabat_minS 70`

#### `--metabat_maxEdges`
- **Type**: Integer
- **Default**: `200`
- **Description**: Maximum edges in the MetaBAT2 graph
- **Example**: `--metabat_maxEdges 250`

#### `--metabat_pTNF`
- **Type**: Integer
- **Default**: `0`
- **Description**: TNF probability threshold for MetaBAT2
- **Example**: `--metabat_pTNF 1`

#### `--metabat_minCV`
- **Type**: Integer
- **Default**: `1`
- **Description**: Minimum coefficient of variation for MetaBAT2
- **Example**: `--metabat_minCV 2`

#### `--metabat_minCVSum`
- **Type**: Integer
- **Default**: `1`
- **Description**: Minimum sum of coefficient of variation for MetaBAT2
- **Example**: `--metabat_minCVSum 2`

#### `--metabat_minClsSize`
- **Type**: Integer
- **Default**: `200000`
- **Description**: Minimum cluster size (bp) for MetaBAT2
- **Example**: `--metabat_minClsSize 250000`

---

## ARG Prediction Parameters

### `--deeparg_min_prob`
- **Type**: Number
- **Default**: `0.8`
- **Range**: 0-1
- **Description**: Minimum probability threshold for DeepARG predictions
- **Example**: `--deeparg_min_prob 0.9`

### `--deeparg_arg_alignment_identity`
- **Type**: Integer
- **Default**: `50`
- **Range**: 0-100
- **Description**: Minimum alignment identity percentage for DeepARG
- **Example**: `--deeparg_arg_alignment_identity 60`

### `--deeparg_arg_alignment_evalue`
- **Type**: String
- **Default**: `1e-10`
- **Description**: Maximum E-value for DeepARG alignments
- **Example**: `--deeparg_arg_alignment_evalue 1e-15`

### `--deeparg_arg_alignment_overlap`
- **Type**: Number
- **Default**: `0.8`
- **Range**: 0-1
- **Description**: Minimum alignment overlap for DeepARG
- **Example**: `--deeparg_arg_alignment_overlap 0.9`

### `--deeparg_arg_num_alignments_per_entry`
- **Type**: Integer
- **Default**: `1000`
- **Description**: Number of alignments per entry for DeepARG
- **Example**: `--deeparg_arg_num_alignments_per_entry 1500`

### `--deeparg_model_version`
- **Type**: String
- **Default**: `v2`
- **Options**: `v1`, `v2`
- **Description**: DeepARG model version
- **Example**: `--deeparg_model_version v2`

### `--rgi_card_version`
- **Type**: String
- **Default**: `latest`
- **Description**: CARD database version for RGI. Use 'latest' or specify a version number
- **Example**: `--rgi_card_version latest` or `--rgi_card_version 3.2.9`

### `--rgi_include_wildcard`
- **Type**: Boolean
- **Default**: `true`
- **Description**: Include WildCARD variants for extended allelic diversity (recommended for environmental samples). The k-mer database of the pathogen-of-origin step (`rgi kmer_query`) comes with WildCARD, so with the downloaded database and `false` that step is skipped (a warning says so)
- **Example**: `--rgi_include_wildcard true`

### `--rgi_aligner`
- **Type**: String
- **Default**: `kma`
- **Options**: `kma`, `bowtie2`, `bwa`
- **Description**: Read aligner for RGI bwt. KMA is recommended for optimal performance with CARD
- **Example**: `--rgi_aligner kma`

### `--rgi_kmer_size`
- **Type**: Integer
- **Default**: `61`
- **Description**: K-mer size for pathogen-of-origin prediction
- **Example**: `--rgi_kmer_size 61`

### `--rgi_min_kmer_coverage`
- **Type**: Integer
- **Default**: `10`
- **Description**: Minimum k-mer coverage threshold for pathogen prediction
- **Example**: `--rgi_min_kmer_coverage 10`

### `--rgi_summary_min_coverage`
- **Type**: Number
- **Default**: `80`
- **Range**: 0-100
- **Description**: `rgi bwt` lists every gene with at least one mapped read; it has no detection rule. In the summaries an RGI gene counts when reads cover at least this % of its best-covered reference sequence (`Percent Coverage` in the allele table): the `*_cov<N>` columns and the plots of `RGI_REPORT`, and the RGI calls in `00_summary/`. The raw RGI tables and the existing `RGI_summary_report.csv` columns are not filtered
- **Example**: `--rgi_summary_min_coverage 60`

### `--kmarvel_min_coverage`
- **Type**: Integer
- **Default**: `90`
- **Range**: 1-100
- **Description**: K-MARVEL minimum coverage (%) of the CARD reference protein by an assembled ARG (`--coverage`)
- **Example**: `--kmarvel_min_coverage 80`

### `--kmarvel_min_identity`
- **Type**: Integer
- **Default**: `90`
- **Range**: 1-100
- **Description**: K-MARVEL minimum identity (%) of an assembled ARG to its CARD reference protein (`--identity`)
- **Example**: `--kmarvel_min_identity 95`

### `--kmarvel_min_read_length`
- **Type**: Integer
- **Default**: `100`
- **Range**: 1-199
- **Description**: K-MARVEL minimum read length in bp (`--min-length`); reads shorter than this after fastp trimming are ignored. K-MARVEL runs in short-read mode, which rejects values of 200 or more; its k-mer (6) and extension (9) sizes are fixed by that mode and not exposed
- **Example**: `--kmarvel_min_read_length 75`

---

## Functional Annotation Parameters

### `--featurecounts_multimap`
- **Type**: String
- **Default**: `primary`
- **Options**: `primary`, `all`, `none`
- **Description**: Multi-mapping policy for featureCounts gene quantification (`--contig_level_functional` branch, published to `07_functional_annotation/gene_abundance/`): `primary` counts primary alignments only, `all` counts every reported alignment (featureCounts `-M`), `none` excludes multi-mapping reads entirely. Counting is read-level (each mate counted separately), not fragment-level
- **Example**: `--featurecounts_multimap all`
- **Note**: Both aligners keep one record per read (bwa-mem3: secondary and supplementary records are dropped; bowtie2 reports one alignment per read by default), so the three settings coincide in practice; the parameter makes the counting policy explicit

### `--microbecensus`
- **Type**: Boolean
- **Default**: `true`
- **Description**: Run MicrobeCensus on the host-removed reads (when `--contig_level_functional` or `--read_level_functional` is enabled; published to `07_functional_annotation/microbecensus/`, once per sample even with both branches on) to estimate average genome size and genome equivalents, enabling copies-per-genome-equivalent (CPGE) normalization in the `07_functional_annotation/summary/` tables of both branches (the SUPER-FOCUS read backend has no CPGE — a SEED hit carries no gene length — so there its AGS/GE are only reported in `read_sample_summary.tsv`; the HUMAnN backend uses them for CPGE = RPK / GE)
- **Example**: `--microbecensus false`
- **Note**: Failure is non-fatal by design: a sample whose MicrobeCensus run fails (reads under 50 bp, too few marker-gene hits, or an estimate outside the 0.5–20 Mb plausibility window) falls back to TPM-only with empty `cpge` fields, recorded as `status = unavailable` in `summary/ags_and_ge.tsv`. Estimates need a few hundred thousand reads to be meaningful

### `--functional_cazy`
- **Type**: Boolean
- **Default**: `true`
- **Options**: `true`, `false`
- **Description**: Run run_dbcan v5 CAZy annotation on the predicted proteins (`--contig_level_functional` branch, published to `07_functional_annotation/dbcan/`): DIAMOND vs CAZy plus pyHMMER vs the dbCAN and dbCAN-sub HMM databases, consolidated into a per-gene `overview.tsv` with the per-tool calls retained, and dbCAN-sub substrate predictions in `dbCANsub_hmm_results.tsv`. The calls feed the summary tables as `db = dbcan` / `backend = run_dbcan` rows alongside the eggNOG-derived CAZy calls, plus dedicated `function_wide_cazy_dbcan_{tpm,cpge}.tsv` matrices (the two CAZy backends are never merged into one matrix — a gene called by both would double-count). Downloads the dbCAN database (~7.4 GB, `--dbcan_db`) on first use
- **Example**: `--functional_cazy false`
- **Note**: With `--functional_cazy false` the summary tables are eggNOG-only: the `cazy_dbcan` wide matrices are still written but header-only, and the `cazy_dbcan` rows in `annotated_fraction.tsv` read zero

### `--dbcan_consensus`
- **Type**: String
- **Default**: `recommended`
- **Options**: `recommended`, `any`
- **Description**: Which dbCAN calls feed the summary tables: `recommended` uses the tool's `Recommend Results` column (calls supported by at least 2 of DIAMOND / dbCAN HMM / dbCAN-sub — the dbCAN authors' guidance), `any` uses the union of the per-tool calls. The published `overview.tsv` always retains all per-tool columns regardless of this setting
- **Example**: `--dbcan_consensus any`

### `--woltka_uniq`
- **Type**: Boolean
- **Default**: `false`
- **Options**: `true`, `false`
- **Description**: Multi-mapping policy of the Woltka read-level backend. Reads are aligned to WoLr2 with the SHOGUN multi-hit Bowtie2 settings (`--very-sensitive -k 16 --np 1 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.05`, the WoL/Qiita standard). By default Woltka divides a read whose reported alignments (up to 16 within the score threshold) overlap k ORFs 1/k to each, so the read totals stay equal to the number of reads assigned. With `true` such ambiguous reads are left unassigned instead (Woltka `--uniq`); how many is reported per sample in `read_sample_summary.tsv` (`reads_unassigned_ambiguous`)
- **Example**: `--woltka_uniq true`

### `--superfocus_aligner`
- **Type**: String
- **Default**: `diamond`
- **Options**: `diamond`, `mmseqs2`
- **Description**: Search backend of the SUPER-FOCUS read-level backend (`--read_level_functional superfocus`): `diamond` (DIAMOND 2.2.1 blastx) or `mmseqs2` (MMseqs2 18 easy-search). Selects which DB_90 archive is downloaded; a `--custom_superfocus_db` must contain that aligner's files. MMseqs2 builds its k-mer index for the database at every run (~13 GB RAM, ~45 s on 16 CPUs) and, in its fast mode, its result for borderline reads can vary slightly with the thread count; DIAMOND results are stable. In a single spot check MMseqs2's fast mode was also somewhat less sensitive (one 2,000-read test sample: DIAMOND 1,011 reads hit, MMseqs2 924). Identity / alignment-length / e-value thresholds are set in the `SUPERFOCUS` `ext.args` (`config/modules.config`), not as parameters
- **Example**: `--superfocus_aligner mmseqs2`

---

## Alignment Parameters

### `--read_aligner`
- **Type**: String
- **Default**: `bwamem3`
- **Options**: `bwamem3`, `bowtie2`
- **Description**: Short-read aligner for host/phiX removal, reads-to-contigs mapping (binning depth, contig coverage, gene counting) and co-assembly bin depth. `bwamem3` uses [bwa-mem3](https://github.com/fg-labs/bwa-mem3) 0.14.0. It replaced `bowtie2` as the default after a benchmark on five real metagenomes ([docs/benchmarks](benchmarks/bwamem3_vs_bowtie2.md)): host removal 4.1× faster, reads→contigs mapping 3.9–6.2× faster with the alignment filter, and the same bins and GTDB species. `bowtie2` stays available as a fallback for CPUs without AVX2 and for nodes without ~75 GB RAM for the human index build. Woltka and HUMAnN always use bowtie2 with their own databases. Differences from `bowtie2`:
  - **Decontamination is strict.** A pair counts as clean only when both mates are unmapped. bowtie2 `--un-conc` kept every pair without a concordant alignment, so a pair with a single host mate survived. The `--bowtie_*` parameters do not apply; `--decontam_min_score` sets the threshold instead.
  - **Decontamination index.** It is built from the same genome FASTA files and published to `<databases_dir>/bwamem3_index/contaminants_index`. Building it for a human host needs **~75 GB RAM** (a few minutes). The request is sized from the genome FASTA files (~24 bytes per base, a gzipped file counted at 4× its size; ~90 GB for the default gzipped human genome), doubles on retry and is capped by `--max_memory`. Mapping against it peaked at ~16 GB per task in the benchmark, against ~4 GB for bowtie2 (the task requests 72 GB, `process_high`, capped by `--max_memory`).
  - **Alignments to contigs.** Secondary and supplementary records are dropped, so every read keeps one record as with bowtie2. Short or divergent local hits are turned into unmapped reads (`--bwamem3_min_identity`, `--bwamem3_min_aligned_length`). The co-assembly contig index is built once and shared by all alignments.
  - **Co-assembly bin depth.** Each sample is mapped once against all refined bins together instead of once per bin.
  - **CPU.** The bioconda bwa-mem3 build needs a CPU with AVX2 or newer (x86-64-v3, Haswell 2013+). On an older CPU the `BWAMEM3_*` tasks stop with "Illegal instruction" (exit status 132): use `--read_aligner bowtie2`.
- **Example**: `--read_aligner bowtie2`

### `--bwamem3_dense_sa`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Build the bwa-mem3 decontamination index with a denser suffix array (`-u 1`, one row in two instead of one in eight). Mapping is faster, at ~12 GB more RAM per decontamination task for a human host. Only for `--read_aligner bwamem3`.
- **Example**: `--bwamem3_dense_sa true`

### `--decontam_min_score`
- **Type**: Integer
- **Default**: `30`
- **Description**: bwa-mem3 minimum alignment score (`-T`) for a read to count as a contaminant hit. Only for `--read_aligner bwamem3`. With bwa's default scoring, 30 corresponds to ~30 matching bases. The threshold applies to each mate on its own: decontamination runs bwa-mem3 with pairing off (`-P`), because in paired mode bwa mem skips `-T` for a chosen pair and random ~20 bp host hits on overlapping mates would remove microbial reads.
- **Example**: `--decontam_min_score 40`

### `--bwamem3_min_identity`
- **Type**: Number (0–1)
- **Default**: `0.95`
- **Description**: Minimum identity, `1 - NM / (M + I + D)`, for a bwa-mem3 read-to-contig or read-to-bin alignment to count. A failing read is kept as an unmapped read. Only for `--read_aligner bwamem3`. It applies to contig depth, binning, featureCounts and co-assembly bin depth, not to decontamination. bwa mem aligns locally, so reads from genomes the assembly lacks still place short or divergent partial hits. A MAPQ cut does not remove them, because such a hit is often the read's only placement. In the bowtie2 vs bwa-mem3 benchmark ([docs/benchmarks](benchmarks/bwamem3_vs_bowtie2.md)), these hits inflated gene counts and bin depth in poorly covered samples. `0` turns the identity check off.
- **Example**: `--bwamem3_min_identity 0.97`

### `--bwamem3_min_aligned_length`
- **Type**: Integer
- **Default**: `50`
- **Description**: Minimum number of aligned bases (CIGAR `M`) for a bwa-mem3 read-to-contig or read-to-bin alignment to count; a failing read is kept as an unmapped read. Only for `--read_aligner bwamem3`; see `--bwamem3_min_identity`. Reads clipped at a contig end still count when 50 bases align. A read shorter than the threshold (a 2×50 library, or reads fastp trimmed below it) needs 80% of its length aligned instead, so short reads that align well still count while a short local hit on a long read does not. `0` turns the length check off.
- **Example**: `--bwamem3_min_aligned_length 40`

Bowtie2 parameters for read alignment during host filtering (`--read_aligner bowtie2`):

### `--bowtie_ma`
- **Type**: Integer
- **Default**: `2`
- **Description**: Match bonus for Bowtie2
- **Example**: `--bowtie_ma 3`

### `--bowtie_mp`
- **Type**: String
- **Default**: `6,2`
- **Description**: Mismatch penalty for Bowtie2 (max, min)
- **Example**: `--bowtie_mp "7,3"`

### `--bowtie_score_min`
- **Type**: String
- **Default**: `G,15,6`
- **Description**: Minimum alignment score function for Bowtie2
- **Example**: `--bowtie_score_min "G,20,8"`

### `--bowtie_k`
- **Type**: Integer
- **Default**: `1`
- **Description**: Number of alignments to report for Bowtie2
- **Example**: `--bowtie_k 2`

### `--bowtie_N`
- **Type**: Integer
- **Default**: `1`
- **Description**: Number of mismatches allowed in seed for Bowtie2
- **Example**: `--bowtie_N 0`

### `--bowtie_L`
- **Type**: Integer
- **Default**: `20`
- **Description**: Seed length for Bowtie2
- **Example**: `--bowtie_L 22`

### `--bowtie_R`
- **Type**: Integer
- **Default**: `2`
- **Description**: Number of re-seeding attempts for Bowtie2
- **Example**: `--bowtie_R 3`

### `--bowtie_i`
- **Type**: String
- **Default**: `S,1,0.75`
- **Description**: Interval function for Bowtie2 seeding
- **Example**: `--bowtie_i "S,1,0.5"`

---

## Resource Limit Parameters

Every task's request is capped at these values. With the local executor, `max_cpus` and
`max_memory` are also the run's machine-wide pool. Set them on the command line or in a
`-params-file` to change both. A `-c` config file that sets them in `params {}` changes the
per-task caps but not the pool, so it must also set `executor { cpus = ...; memory = ... }`.
Otherwise a task asking for more than the pool fails with `Process requirement exceeds
available CPUs` (or `memory`); see manual §10.1. `-profile test` sets 2 CPUs / 6 GB / 6 h.

### `--max_cpus`
- **Type**: Integer
- **Default**: `16`
- **Description**: Maximum number of CPUs that can be requested for any single job; also the local executor's CPU pool (and queue size)
- **Example**: `--max_cpus 32`

### `--max_memory`
- **Type**: String
- **Default**: `128.GB`
- **Description**: Maximum amount of memory that can be requested for any single job; also the local executor's memory pool
- **Example**: `--max_memory 256.GB`

### `--max_time`
- **Type**: String
- **Default**: `240.h`
- **Description**: Maximum amount of time that can be requested for any single job
- **Example**: `--max_time 480.h`

---

## Advanced Parameters

### `--help`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Display help text and exit
- **Example**: `--help`

### `--version`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Display version and exit
- **Example**: `--version`

> **Note**: Parameters are validated against `nextflow_schema.json` at startup. Any parameter not declared by the pipeline — including misspelled or removed ones — aborts the run with an "unrecognised parameter" error.

---

## Cloud & HPC Profile Parameters

These parameters are read by the execution profiles in `conf/` (`-profile aws`, `gcp`, `azure`, `slurm`) and are ignored elsewhere. The `--*_workdir` parameter of the chosen cloud profile is **required** — the pipeline exits immediately with an error if it is missing.

### AWS Batch (`-profile aws`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--aws_workdir` | *(required)* | S3 work directory, e.g. `s3://my-bucket/work` |
| `--aws_queue` | `default` | AWS Batch job queue |
| `--aws_region` | `us-east-1` | AWS region |
| `--aws_cli_path` | `/home/ec2-user/miniconda/bin/aws` | Path to the `aws` CLI on the Batch AMI |
| `--aws_job_role` | — | IAM job role ARN |
| `--aws_execution_role` | — | IAM execution role ARN |
| `--aws_volumes` | — | Container volumes for scratch space |

### Google Cloud Batch (`-profile gcp`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gcp_workdir` | *(required)* | GCS work directory, e.g. `gs://my-bucket/work` |
| `--gcp_project` | — | GCP project id |
| `--gcp_region` | `us-central1` | GCP location |
| `--gcp_spot` | `false` | Use spot/preemptible instances |
| `--gcp_boot_disk_size` | `50.GB` | Boot disk size per VM |
| `--gcp_service_account` | — | Service account email |
| `--gcp_network` / `--gcp_subnetwork` | — | VPC network settings |
| `--gcp_private_address` | `false` | Use private IP addresses only |

### Azure Batch (`-profile azure`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--azure_workdir` | *(required)* | Blob work directory, e.g. `az://my-container/work` |
| `--azure_pool` | `auto` | Batch pool name |
| `--azure_region` | `eastus` | Azure location |
| `--azure_batch_account` / `--azure_batch_key` | — | Batch account credentials |
| `--azure_storage_account` / `--azure_storage_key` | — | Storage account credentials |
| `--azure_delete_pools` | `true` | Delete auto pools on completion |
| `--azure_vm_type` | `Standard_D4_v3` | VM type for auto pools |
| `--azure_vm_count` | `1` | Initial VM count |
| `--azure_max_vms` | `10` | Maximum VM count |
| `--azure_spot` | `false` | Use low-priority (spot) VMs |

### SLURM (`-profile slurm`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--slurm_queue` | `normal` | Default partition |
| `--slurm_queue_short` / `--slurm_queue_long` | — | Partitions for short/long jobs |
| `--slurm_account` | — | Account for job submission |
| `--slurm_queue_size` | `100` | Executor queue size |
| `--singularity_cache` | — | Singularity image cache directory |
| `--singularity_run_options` | — | Extra singularity run options |

---

## Parameter Usage Examples

### Minimal Run
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker
```

### Full Pipeline with Custom Parameters
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --quality_control true \
    --assembly_mode assembly \
    --taxonomic_profiler kraken2 \
    --include_binning true \
    --binners semibin,metabat2 \
    --read_arg_prediction true \
    --contig_tax_and_arg true \
    --metawrap_completeness 70 \
    --metawrap_contamination 5 \
    --semibin_env_model human_gut \
    --max_cpus 32 \
    --max_memory 256.GB \
    -profile singularity
```

### Using Custom Databases
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --custom_kraken_db /shared/db/kraken2 \
    --custom_checkm2_db /shared/db/checkm2.dmnd \
    --custom_gtdbtk_db /shared/db/gtdbtk/release232 \
    --databases_dir /shared/databases \
    -profile docker
```

---

*BugBuster v2.0.0dev - Complete Parameter Reference*
