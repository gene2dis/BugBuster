# Changelog

All notable changes to BugBuster will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
  - Added `versions.yml` output for software tracking
  - Added `stub` blocks for dry-run testing
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
