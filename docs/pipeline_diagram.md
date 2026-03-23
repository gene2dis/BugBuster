# BugBuster Pipeline Diagram

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              BugBuster Pipeline                                  │
└─────────────────────────────────────────────────────────────────────────────────┘

                                 ┌─────────────┐
                                 │  Input CSV  │
                                 │ Samplesheet │
                                 └──────┬──────┘
                                        │
                    ┌───────────────────┼───────────────────┐
                    │                   │                   │
                    ▼                   ▼                   ▼
          ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐
          │  INPUT_CHECK    │  │  PREPARE_       │  │                 │
          │  (Validation)   │  │  DATABASES      │  │                 │
          └────────┬────────┘  └────────┬────────┘  │                 │
                   │                    │            │                 │
                   └────────────────────┼────────────┘                 │
                                        │                              │
                                        ▼                              │
                              ┌─────────────────┐                      │
                              │   QC SUBWF      │◄─────────────────────┘
                              │  (FASTP, Host/  │  (Databases: PhiX,
                              │   PhiX Filter)  │   Host indexes)
                              └────────┬────────┘
                                       │
                                       ▼
                             ┌─────────────────┐
                             │  Clean Reads    │
                             └────────┬────────┘
                                      │
        ┌─────────────────────────────┼─────────────────────────────┐
        │                             │                             │
        ▼                             ▼                             ▼
┌─────────────────┐         ┌─────────────────┐         ┌─────────────────┐
│   TAXONOMY      │         │   ASSEMBLY      │         │  READ ARG       │
│   SUBWF         │         │   SUBWF         │         │  PREDICTION     │
│ (if enabled)    │         │ (if enabled)    │         │ (if enabled)    │
└────────┬────────┘         └────────┬────────┘         └────────┬────────┘
         │                           │                           │
         ▼                           ▼                           ▼
┌─────────────────┐         ┌─────────────────┐         ┌─────────────────┐
│ KRAKEN2/BRACKEN │         │    MEGAHIT      │         │  KARGVA→KARGA   │
│   SOURMASH      │         │     BBMAP       │         │   ARGS-OAP      │
│   PHYLOSEQ      │         │  BOWTIE2/SAM    │         │  ARG_NORM_RPT   │
└─────────────────┘         └────────┬────────┘         └─────────────────┘
                                     │
                 ┌───────────────────┼───────────────────┐
                 │                   │                   │
                 ▼                   ▼                   ▼
        ┌───────────────┐   ┌───────────────┐   ┌───────────────┐
        │   BINNING     │   │ CONTIG TAX &  │   │ METACERBERUS  │
        │   SUBWF       │   │  ARG ANALYSIS │   │  (assembly    │
        │ (if enabled)  │   │ (if enabled)  │   │   mode only)  │
        └───────┬───────┘   └───────┬───────┘   └───────────────┘
                │                   │
                ▼                   ▼
        ┌───────────────┐        ┌───────────────┐
        │ Selected      │   │  NT_BLASTN    │
        │ Binners:      │   │  SAMTOOLS_IDX │
        │ METABAT2 *    │   │  BLOBTOOLS    │
        │ SEMIBIN  *    │   │  BLOBPLOT     │
        │ COMEBIN  *    │   └───────┬───────┘
        │ (* if chosen) │           │
        └───────┬───────┘           ▼
                │           ┌───────────────┐
                ▼           │  PRODIGAL     │
        ┌───────────────┐   │  DEEPARG      │
        │   METAWRAP    │   │  ARG_CONTIG   │
        │(if ≥2 binners)│   │  ARG_BLOBPLOT │
        └───────┬───────┘   └───────────────┘
                │
        ┌───────┴───────┐
        │               │
        ▼               ▼
┌───────────────┐ ┌───────────────┐
│   CHECKM2     │ │   GTDB-TK     │
│  (Quality)    │ │  (Taxonomy)   │
└───────┬───────┘ └───────┬───────┘
        │               │
        └───────┬───────┘
                ▼
        ┌───────────────┐
        │  BIN REPORTS  │
        │ (Quality/Tax/ │
        │   Summary)    │
        └───────────────┘
                │
                ▼
        ┌───────────────┐
        │   MULTIQC     │
        │  (Aggregate)  │
        └───────────────┘
```

## Subworkflow & Module Descriptions

### Subworkflows
| Subworkflow | Description |
|-------------|-------------|
| `INPUT_CHECK` | Validate samplesheet format and file existence |
| `PREPARE_DATABASES` | Download and format all required databases |
| `QC` | Quality control and host decontamination |
| `TAXONOMY` | Taxonomic profiling with Kraken2 or Sourmash |
| `ASSEMBLY` | Metagenome assembly (per-sample or co-assembly) |
| `BINNING` | Unified binning workflow (mode-agnostic) |

### Quality Control (QC Subworkflow)
| Module | Description |
|--------|-------------|
| `FASTP` | Read quality filtering and adapter trimming |
| `QFILTER` | Extract QC reports and filter by read count |
| `BOWTIE2` (Host) | Remove host genome reads |
| `BOWTIE2` (PhiX) | Remove PhiX spike-in contamination |
| `READS_REPORT` | Generate read count summary report |

### Taxonomic Profiling
| Module | Description |
|--------|-------------|
| `KRAKEN2` | K-mer based taxonomic classification |
| `BRACKEN` | Abundance estimation from Kraken2 |
| `SOURMASH` | MinHash-based taxonomic profiling |
| `PHYLOSEQ` | Generate R Phyloseq objects |

### Assembly
| Module | Description |
|--------|-------------|
| `MEGAHIT` | De novo metagenome assembly |
| `BBMAP` | Contig length filtering |

### Binning (BINNING Subworkflow - Mode-Agnostic)
| Module | Condition | Description |
|--------|-----------|-------------|
| `CALCULATE_DEPTH` | if `metabat2` selected | Calculate contig depth from BAM files |
| `METABAT2` | if `metabat2` in `--binners` | Binning by coverage and composition |
| `SEMIBIN` | if `semibin` in `--binners` | Semi-supervised binning |
| `COMEBIN` | if `comebin` in `--binners` | Contrastive learning binning |
| `METAWRAP` | if ≥2 binners selected | Bin refinement and consolidation |
| `BOWTIE2_SAMTOOLS_DEPTH` | co-assembly only | Calculate bin coverage |
| `BEDTOOLS` | co-assembly only | Coverage statistics |

### Quality Assessment & Reporting
| Module | Description |
|--------|-------------|
| `CHECKM2` | Bin completeness and contamination |
| `GTDB-TK` | Bin taxonomic classification |
| `BIN_QUALITY_REPORT` | Generate bin quality report (assembly mode) |
| `BIN_TAX_REPORT` | Generate bin taxonomy report (assembly mode) |
| `BIN_SUMMARY` | Comprehensive bin summary (co-assembly mode) |

### ARG Prediction
| Module | Description |
|--------|-------------|
| `KARGA` | Read-level ARG detection |
| `KARGVA` | Read-level ARG variant detection |
| `DEEPARG` | Contig-level ARG prediction |
| `ARGS-OAP` | ARG normalization |

### Contig Analysis
| Module | Description |
|--------|-------------|
| `NT_BLASTN` | Contig taxonomic assignment via BLAST |
| `SAMTOOLS_INDEX` | Index BAM files for BLOBTOOLS |
| `BLOBTOOLS` | Generate blob tables for contig taxonomy |
| `BLOBPLOT` | Visualize contig taxonomy |
| `PRODIGAL_CONTIGS` | ORF prediction on contigs |
| `DEEPARG_CONTIGS` | ARG prediction on contigs |
| `ARG_CONTIG_LEVEL_REPORT` | Generate contig-level ARG report |
| `ARG_BLOBPLOT` | Visualize ARG distribution on contigs |
| `METACERBERUS` | Functional annotation (assembly mode only) |

### Bin ARG Analysis
| Module | Description |
|--------|-------------|
| `PRODIGAL_BINS` | ORF prediction on bins |
| `DEEPARG_BINS` | ARG prediction on bins |
| `ARG_FASTA_FORMATTER` | Format ARG sequences |
| `CLUSTERING` | Cluster ARG sequences |

## Execution Modes

### Assembly Mode (`--assembly_mode`)

| Mode | Description |
|------|-------------|
| `assembly` | Each sample assembled independently |
| `coassembly` | All samples pooled and assembled together |
| `none` | Skip assembly, taxonomy only |

### Taxonomic Profiler (`--taxonomic_profiler`)

| Option | Description |
|--------|-------------|
| `kraken2` | Use Kraken2 + Bracken |
| `sourmash` | Use Sourmash gather |
| `none` | Skip taxonomic profiling |

## Key Pipeline Features

### Unified Binning Workflow
The BINNING subworkflow is **mode-agnostic**, handling both per-sample assembly and co-assembly modes with a single unified workflow:
- Takes `contigs_and_bam` input directly
- Runs only the binners specified in `--binners` (default: `semibin`)
- **MetaWRAP is only invoked when ≥2 binners are selected** — when a single binner is used, its output passes directly downstream
- Assesses quality with CheckM2 and classifies with GTDB-Tk (always runs)
- Generates mode-specific reports:
  - **Assembly mode**: Simple quality and taxonomy reports
  - **Co-assembly mode**: Additional bin coverage analysis and comprehensive summary

**Binner selection examples:**
```bash
# Default - single fast binner, no MetaWRAP
--include_binning true --binners semibin

# Two binners - triggers MetaWRAP refinement
--include_binning true --binners semibin,metabat2

# All three binners
--include_binning true --binners semibin,metabat2,comebin
```

### Database Preparation
The PREPARE_DATABASES subworkflow handles all database downloads and formatting:
- Kraken2/Sourmash databases for taxonomy
- PhiX and host indexes for QC
- DeepARG, BLAST, and taxdump for contig analysis
- CheckM2 and GTDB-Tk databases for binning
- KARGA/KARGVA databases for read-level ARG prediction
