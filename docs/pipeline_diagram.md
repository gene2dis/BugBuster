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
                                        ▼
                              ┌─────────────────┐
                              │  INPUT_CHECK    │
                              │  (Validation)   │
                              └────────┬────────┘
                                       │
                    ┌──────────────────┼──────────────────┐
                    │                  │                  │
                    ▼                  ▼                  ▼
           ┌───────────────┐  ┌───────────────┐  ┌───────────────┐
           │   FASTP       │  │   BOWTIE2     │  │   BOWTIE2     │
           │ (QC Filter)   │  │ (PhiX Remove) │  │ (Host Remove) │
           └───────┬───────┘  └───────┬───────┘  └───────┬───────┘
                   │                  │                  │
                   └──────────────────┼──────────────────┘
                                      │
                                      ▼
                            ┌─────────────────┐
                            │  Clean Reads    │
                            └────────┬────────┘
                                     │
           ┌─────────────────────────┼─────────────────────────┐
           │                         │                         │
           ▼                         ▼                         ▼
  ┌─────────────────┐      ┌─────────────────┐      ┌─────────────────┐
  │    TAXONOMY     │      │    ASSEMBLY     │      │  ARG PREDICTION │
  │  (if enabled)   │      │  (if enabled)   │      │  (if enabled)   │
  └────────┬────────┘      └────────┬────────┘      └────────┬────────┘
           │                        │                        │
           ▼                        ▼                        ▼
  ┌─────────────────┐      ┌─────────────────┐      ┌─────────────────┐
  │ KRAKEN2/SOURMASH│      │    MEGAHIT      │      │  KARGA/KARGVA   │
  │    BRACKEN      │      │     BBMAP       │      │   ARGS-OAP      │
  │   PHYLOSEQ      │      │  (Filter)       │      └─────────────────┘
  └─────────────────┘      └────────┬────────┘
                                    │
                    ┌───────────────┼───────────────┐
                    │               │               │
                    ▼               ▼               ▼
           ┌───────────────┐ ┌───────────────┐ ┌───────────────┐
           │   BINNING     │ │ CONTIG TAX    │ │ METACERBERUS  │
           │ (if enabled)  │ │ (if enabled)  │ │ (if enabled)  │
           └───────┬───────┘ └───────┬───────┘ └───────────────┘
                   │                 │
                   ▼                 ▼
           ┌───────────────┐ ┌───────────────┐
           │  METABAT2     │ │   BLASTN      │
           │  SEMIBIN      │ │  BLOBTOOLS    │
           │  COMEBIN      │ │   DEEPARG     │
           └───────┬───────┘ └───────────────┘
                   │
                   ▼
           ┌───────────────┐
           │   METAWRAP    │
           │ (Refinement)  │
           └───────┬───────┘
                   │
           ┌───────┴───────┐
           │               │
           ▼               ▼
   ┌───────────────┐ ┌───────────────┐
   │   CHECKM2     │ │   GTDB-TK     │
   │  (Quality)    │ │  (Taxonomy)   │
   └───────────────┘ └───────────────┘
                   │
                   ▼
           ┌───────────────┐
           │    REPORTS    │
           └───────────────┘
```

## Module Descriptions

### Input Processing
| Module | Description |
|--------|-------------|
| `INPUT_CHECK` | Validate samplesheet format and file existence |

### Quality Control (QC)
| Module | Description |
|--------|-------------|
| `FASTP` | Read quality filtering and adapter trimming |
| `BOWTIE2` (PhiX) | Remove PhiX spike-in contamination |
| `BOWTIE2` (Host) | Remove host genome reads |

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

### Binning
| Module | Description |
|--------|-------------|
| `METABAT2` | Binning by coverage and composition |
| `SEMIBIN` | Semi-supervised binning |
| `COMEBIN` | Contrastive learning binning |
| `METAWRAP` | Bin refinement and consolidation |

### Quality Assessment
| Module | Description |
|--------|-------------|
| `CHECKM2` | Bin completeness and contamination |
| `GTDB-TK` | Bin taxonomic classification |

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
| `BLASTN` | Contig taxonomic assignment |
| `BLOBTOOLS` | Contig visualization |
| `PRODIGAL` | ORF prediction |
| `METACERBERUS` | Functional annotation |

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
