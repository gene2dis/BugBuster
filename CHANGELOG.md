# Changelog

All notable changes to BugBuster will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

The pipeline self-reports this state as `1.1.0dev` (manifest version) until the
next release is tagged.

### Added

- **MAG-level functional annotation (`--mag_level_functional`)**
  - Bakta 1.12.1 annotation of every MetaWRAP-refined bin, one task per bin
    (nf-core `bakta/bakta` module, patched to emit `versions.yml`; per-bin
    fan-out of the refined-bins directory), published to
    `07_functional_annotation/mags/<sample>/` as `<sample>_<bin>.{gff3,gbff,faa,fna,tsv,txt,hypotheticals.tsv,hypotheticals.faa}`
  - Requires `--include_binning` and at least two `--binners` (validated at
    launch): the MetaWRAP completeness/contamination filter only runs with
    ≥2 binners, and only quality-filtered bins are annotated. Bins are
    expected to be bacterial; the pipeline does not exclude
    archaeal/eukaryotic/viral bins (documented)
  - Bakta database v6.0 (schema 6) auto-downloaded from the pinned Zenodo
    release: `--bakta_db v6.0-full` (default, 31.9 GB download) or
    `v6.0-light` (1.3 GB; explicit choice, recorded in provenance via the
    `bakta_db` versions entry) / `--custom_bakta_db`. Provisioning refreshes
    the bundled AMRFinderPlus database with `amrfinder_update` (~300 MB; the
    official tarball's copy is too old for the AMRFinderPlus in the pinned
    bakta container and would fail every annotation at the AMR expert step —
    custom databases need the same one-time refresh, see troubleshooting)
  - Independent of `--contig_level_functional` and runs on any container
    engine (docker included) — the singularity/apptainer requirement applies
    only to the contig branch
- **Contig-level functional annotation (`--contig_level_functional`)**
  - Shared Pyrodigal gene calling on contigs (one pass feeds both DeepARG and
    functional annotation), published to `07_functional_annotation/gene_calling/`
  - Two-stage eggNOG-mapper v3 annotation (DIAMOND search + orthology-transfer
    annotation) against the eggNOG 7 database (~44 GB, auto-downloaded;
    `--eggnog_db` / `--custom_eggnog_db`), published to
    `07_functional_annotation/eggnog/`
  - Per-gene abundance quantification with featureCounts over the Pyrodigal
    gene coordinates, with gene ids matching the annotated protein ids;
    per-sample counts in both assembly modes (under co-assembly, dedicated
    per-sample alignments against the co-assembly are added since the pooled
    binning BAM cannot yield per-sample counts), published to
    `07_functional_annotation/gene_abundance/`; multi-mapping policy via
    `--featurecounts_multimap` (`primary`/`all`/`none`)
  - Study-level aggregation (`bin/aggregate_functions.py`) into canonical
    tables published to `07_functional_annotation/summary/`: per-gene
    abundance with TPM and copies per genome equivalent (CPGE), long-format
    gene annotations, per-ontology function abundance (KO, COG, EC, Pfam,
    CAZy; intentional double-counting of multi-term genes), wide TPM and
    CPGE matrices per ontology, a per-sample annotated-fraction report and
    an AGS summary (`ags_and_ge.tsv`); the eggNOG annotations parser is
    version-aware and fails loudly on layout drift
  - MicrobeCensus average genome size estimation on the host-removed reads
    (`--microbecensus`, on by default with the branch), published to
    `07_functional_annotation/microbecensus/`, enabling the CPGE
    normalization. Failure is non-fatal by design: affected samples fall
    back to TPM-only with empty `cpge` fields and `status = unavailable` in
    `ags_and_ge.tsv`. Both published biocontainers of the unmaintained
    upstream tool are broken (Python-3 build: unreleased upstream str/bytes
    fix; Python-2 build: missing libstdc++ for the bundled RAPsearch2); the
    pipeline pins the Python-3 image and routes the call through the
    `bin/run_microbe_census_py3fix.py` shim that patches the one broken
    function
  - run_dbcan v5 CAZy annotation of the predicted proteins
    (`--functional_cazy`, on by default with the branch): protein-mode
    `CAZyme_annotation` (DIAMOND vs CAZy + pyHMMER vs dbCAN/dbCAN-sub HMMs)
    against the pinned dbCAN database release (~7.4 GB, auto-downloaded;
    `--dbcan_db` / `--custom_dbcan_db`), published to
    `07_functional_annotation/dbcan/` with the per-tool overview columns and
    the dbCAN-sub substrate predictions retained. The calls feed the summary
    tables as `db = dbcan` / `backend = run_dbcan` rows alongside the
    eggNOG-derived CAZy calls — reported separately, never merged — plus
    dedicated `function_wide_cazy_dbcan_{tpm,cpge}.tsv` matrices and
    `cazy_dbcan` annotated-fraction rows; the consensus policy is a
    documented parameter (`--dbcan_consensus`: `recommended` = calls
    supported by >= 2 tools, `any` = per-tool union), and the overview
    parser is version-aware like the eggNOG one
  - Note: the branch requires a singularity/apptainer container engine while
    eggNOG-mapper v3 is in beta (no docker image exists upstream); real runs
    under docker/podman abort at launch with an explanatory error
  - `docs/WORKFLOW_DIAGRAM.md` refreshed to cover the functional annotation
    branch (shared Pyrodigal gene calling, eggNOG/dbCAN database preparation,
    per-sample co-assembly counting alignments, the FUNCTIONAL_ANNOTATION
    subworkflow and its parameters)

- **Provenance tracking (audit #27)**
  - `versions.yml` emitted by every live module (previously ~20 in-use modules
    emitted none) and, for the first time, aggregated: each run now writes a
    deduplicated `pipeline_info/software_versions.yml` covering every executed
    process plus the pipeline and Nextflow versions
  - DeepARG database downloads record the tool version and download date; the
    CARD version captured by RGI_LOAD now actually reaches the aggregate report

- **RGI AMR Prediction**
  - New `--rgi_prediction` parameter to enable AMR gene prediction with pathogen-of-origin analysis
  - RGI_LOAD module for automatic CARD database download and preparation
  - RGI_LOAD_WILDCARD module for combining separate CARD and WildCARD databases
  - RGI_BWT module for read-level AMR gene alignment using KMA aligner
  - RGI_KMER module for k-mer based pathogen-of-origin prediction
  - RGI_REPORT module for aggregated results and visualizations
  - Support for WildCARD variants for extended allelic diversity
  - Custom database support via `--custom_rgi_card_db` and `--custom_rgi_wildcard` parameters
  - Flexible database options: automatic download, pre-prepared complete database, or separate CARD + WildCARD
  - Comprehensive documentation in `docs/RGI_WILDCARD_USAGE.md` (originally also `docs/RGI_IMPLEMENTATION_{PLAN,SUMMARY}.md`, since consolidated)
  - RGI-specific parameters: `rgi_card_version`, `rgi_include_wildcard`, `rgi_aligner`, `rgi_kmer_size`, `rgi_min_kmer_coverage`
  - Integration with PREPARE_DATABASES subworkflow for automatic database management
  - Output includes per-sample results, pathogen predictions, and multi-sample summary reports with plots

### Changed

- **Breaking (GTDB-Tk databases)**: GTDB-Tk re-pinned 2.5.2 → 2.7.2 and the
  reference-data registry moved from GTDB R220 to **R232** (~61 GB download,
  pinned release URL; default `--gtdbtk_db release_232`). GTDB-Tk 2.7.x
  accepts only R232 data, so a `databases/gtdbtk/` directory cached by earlier
  pipeline versions (R220) no longer works — delete it to re-download, or pass
  an R232 directory with `--custom_gtdbtk_db`. The upstream-removed
  `--skip_ani_screen` flag was dropped from the `classify_wf` invocation (the
  ANI pre-screen now always runs, against the skani DB bundled in the R232
  package), and the module's synthetic empty-report headers were updated to
  the 2.7 column schema (`fastani_*` → `closest_genome_*`;
  `bin_tax_report.py` already handled both namings)
- **Breaking (output layout)**: contig gene calling switched from the nf-core
  Prodigal module to a shared Pyrodigal step; its outputs moved from
  `05_arg_prediction/contig_level/prodigal/` to
  `07_functional_annotation/gene_calling/` (same predictions, gzipped, now
  produced once for both DeepARG and functional annotation)
- Updated README.md with RGI feature description and usage examples
- Updated docs/manual.md with RGI parameters, output structure, and usage examples
- Updated docs/parameters.md with complete RGI parameter reference
- Enhanced ARG prediction capabilities with complementary tool (RGI alongside KARGA/KARGVA)
- Documentation corrected to match real behavior (audit #26): sample read-count
  filter default, database storage location, output trees, profile lists;
  MultiQC removed from docs, dead wiring, and vendored modules (audit #28;
  re-enabling requires re-vendoring via `nf-core modules install multiqc`)
- Minimum Nextflow version raised from 23.04.0 to **24.04.0** with parameter
  validation now enforced by the nf-schema 2.4.2 plugin (audit #25): unknown
  `--params` are a hard startup error
- Dead documented knobs fixed (audit #14): `fastp_qualified_quality_phred`
  wired, METABAT2 selectors collapsed (pTNF/minCV/minCVSum now delivered),
  `bbmap_lenght` doc typo corrected, unused `mmseqs_*` params removed

### Removed

**Breaking**: parameter validation is now strict (nf-schema
`failUnrecognisedParams`), so passing any removed parameter aborts the run at
startup instead of being silently ignored. Update existing command lines and
`-params-file` YAMLs accordingly.

- `--custom_phiX_index` and `--custom_bowtie_host_index` — replaced by
  `--custom_decontamination_index` (pre-built combined index),
  `--custom_phiX_fasta`, and `--custom_host_fasta` (the pipeline builds the
  combined Bowtie2 index from FASTAs); the old params had been silent no-ops
- `--enable_work_cleanup` — was never consumed; work-dir cleanup is Nextflow's
  `cleanup = true`, enabled by `-profile low_disk`
- `--kraken_db_used`, `--sourmash_db_name` — report database names are derived
  from the selected database
- `--store_filtered_contigs`, `--store_refined_bins` — were no-ops; filtered
  contigs and refined bins are always published
- `--tracedir` — trace/report/timeline/DAG always go to
  `<output>/pipeline_info/`
- `--validationShowHiddenParams`, `--validationSchemaIgnoreParams` —
  nf-validation 1.x options superseded by the nf-schema `validation {}` scope
- All eleven `--mmseqs_*` parameters — unused; clustering settings are fixed
  tiers in `modules/local/clustering`
- Renamed: `--bbmap_lenght` → `--bbmap_length`

### Fixed

- **Contig-branch Pfam terms were fragmented by domain coordinates**: eggNOG-mapper
  v3 writes `PFAMs` values as `<pfam_name>_<start>_<end>`, and the aggregation kept
  each string as the accession, so one Pfam family appeared as many terms (and a gene
  with a repeated domain counted toward several). `summary/` tables now carry the
  plain Pfam name, counted once per gene; an unexpected PFAMs format fails loudly
  (format verified on all 73.9 M pfam values of the eggNOG 7 database)
- Extensive audit-fix series on branch `fix-pending-issues` (2026-08): host
  decontamination DB and `--local` scoring, database download containers and
  script hardening, kraken2/bracken report parsers, contig tax/ARG arm wiring,
  singleton read handling through QC, storeDir misuse, report failure masking,
  silent no-op feature combos, edge-data crashes, publishing gaps (the
  itemized 30-point list and fixing commits are recorded in the per-item
  `fix:`/`feat:` commit messages on that branch)

## [1.0.0] - 2024-01-10

### Added

- **Pipeline Infrastructure**
  - Manifest block with pipeline metadata and Nextflow version requirement (>=23.04.0)
  - Parameter validation schema (`nextflow_schema.json`)
  - Samplesheet validation schema (`assets/schema_input.json`)
  - Input validation subworkflow with file existence checks
  - Dynamic resource allocation with `check_max()` function
  - Execution reports: timeline, trace, and DAG

- **Modular Architecture**
  - Subworkflows for logical pipeline components:
    - `INPUT_CHECK` - Samplesheet validation
    - `PREPARE_DATABASES` - Database download and formatting
    - `QC` - Quality control and host filtering
    - `TAXONOMY` - Taxonomic profiling
    - `ASSEMBLY` - Genome assembly
    - `BINNING` - Metagenomic binning and refinement

- **Cloud Support**
  - AWS Batch configuration with S3 and Fusion filesystem
  - Google Cloud Batch/Life Sciences configuration
  - Azure Batch configuration with auto-scaling pools
  - SLURM HPC configuration with Singularity

- **Container Profiles**
  - Docker profile
  - Singularity profile
  - Podman profile
  - Apptainer profile

- **Module Improvements**
  - Fixed bash null checks using Groovy conditionals
  - Added `versions.yml` output to some modules (full coverage and aggregation
    landed later, in [Unreleased])
  - Added `stub` blocks for some modules (full coverage landed later)
  - Added `meta.yml` descriptors for key modules

- **Documentation**
  - Comprehensive deployment guide (`docs/deployment.md`)
  - Updated README with profiles and cloud examples
  - Contributing guidelines (`CONTRIBUTING.md`)
  - Example samplesheet

### Changed

- Refactored main.nf with improved help/version handling
- Consolidated params block in nextflow.config
- Updated process labels with memory and time definitions
- Improved error handling and retry strategies

### Fixed

- Bash null string comparison issues in FASTP, MEGAHIT, BOWTIE2 modules
- Resource allocation respects max limits on retries

## [0.1.0] - Initial Release

### Added

- Initial pipeline implementation
- Basic QC, assembly, binning workflow
- Kraken2 and Sourmash taxonomic profiling
- ARG prediction with DeepARG, KARGA, KARGVA
- Binning with MetaBAT2, SemiBin, COMEBin
- Bin refinement with MetaWRAP
- Quality assessment with CheckM2
- Taxonomic classification with GTDB-TK

[Unreleased]: https://github.com/gene2dis/BugBuster/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/gene2dis/BugBuster/releases/tag/v1.0.0
[0.1.0]: https://github.com/gene2dis/BugBuster/tree/v0.1.0
