# RGI AMR Gene Prediction Implementation Plan

## Overview

This document outlines the detailed implementation plan for integrating RGI (Resistance Gene Identifier) into the BugBuster pipeline for AMR gene prediction from metagenomic reads with pathogen-of-origin analysis.

## Background

### RGI Tool Capabilities

**RGI** is a tool from the Comprehensive Antibiotic Resistance Database (CARD) that predicts resistomes from protein or nucleotide data, including metagenomics data. Key features:

- **RGI bwt**: Aligns metagenomic reads to CARD AMR alleles using KMA (recommended), Bowtie2, or BWA
- **RGI kmer_query**: K-mer based pathogen-of-origin prediction for AMR genes
- **CARD Database**: Curated AMR reference sequences with experimental validation
- **WildCARD**: Extended allelic variants from in silico predictions across hundreds of pathogens

### Analysis Workflow

1. **Read Alignment** (RGI bwt):
   - Align FASTQ reads to CARD protein homolog models
   - Use KMA aligner (best for redundant databases like CARD)
   - Optional: Include WildCARD variants for broader allelic coverage
   - Output: Allele-level and gene-level mapping results

2. **Pathogen-of-Origin Prediction** (RGI kmer_query):
   - K-mer classification of AMR-positive reads
   - Predict source pathogen for detected AMR genes
   - Requires RGI bwt BAM output as input

## Implementation Architecture

### Module Structure

Following BugBuster's existing patterns (e.g., `@/home/jugalde/pipelines/BugBuster/modules/local/karga/main.nf`, `@/home/jugalde/pipelines/BugBuster/modules/local/deeparg/main.nf`), we will create:

```
modules/local/
├── rgi_bwt/
│   └── main.nf          # RGI read alignment module
├── rgi_kmer/
│   └── main.nf          # Pathogen-of-origin prediction module
└── rgi_load/
    └── main.nf          # Database preparation module
```

### Database Preparation

**Module**: `modules/local/rgi_load/main.nf`

**Purpose**: Download and prepare CARD databases for RGI analysis

**Inputs**:
- Database version parameters (from config)
- Include WildCARD flag

**Process Steps**:
1. Download CARD JSON database
2. Run `rgi card_annotation` to create reference FASTA
3. Load CARD data with `rgi load`
4. Optional: Download and process WildCARD variants
5. Optional: Run `rgi wildcard_annotation` and load
6. Build KMA indices for read alignment

**Outputs**:
- `card_database/`: CARD reference database directory
- `wildcard_database/`: WildCARD variants (optional)
- `kma_index/`: KMA alignment indices

**Configuration Parameters**:
```groovy
params {
    // RGI database options
    rgi_include_wildcard   = true      // Include WildCARD variants
    rgi_card_version       = 'latest'  // CARD version
    rgi_kmer_size          = 61        // K-mer size for pathogen prediction
    
    // Custom database paths (override auto-download)
    custom_rgi_card_json   = null
    custom_rgi_wildcard    = null
}
```

### RGI BWT Module (Read Alignment)

**Module**: `modules/local/rgi_bwt/main.nf`

**Purpose**: Align metagenomic reads to CARD AMR alleles

**Inputs**:
- `tuple val(meta), path(reads)` - Clean paired-end reads from QC
- `path(card_db)` - CARD database directory
- `path(wildcard_db)` - WildCARD database (optional)

**Process Configuration**:
- **Label**: `process_medium` (8 CPUs, 36GB RAM)
- **Container**: `quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0`
- **Aligner**: KMA (default, best performance)

**Command Structure**:
```bash
rgi bwt \
    --read_one ${reads[0]} \
    --read_two ${reads[1]} \
    --aligner kma \
    --threads ${task.cpus} \
    --output_file ${prefix}_rgi_bwt \
    --local \
    ${include_wildcard ? '--include_wildcard' : ''}
```

**Outputs**:
- `tuple val(meta), path("*.allele_mapping_data.txt")` - Allele-level results
- `tuple val(meta), path("*.gene_mapping_data.txt")` - Gene-level results
- `tuple val(meta), path("*.sorted.length_100.bam")` - BAM file for k-mer analysis
- `path("*.overall_mapping_stats.txt")` - Mapping statistics
- `path("*.reference_mapping_stats.txt")` - Reference statistics

**Output Files**:
1. `{sample}_rgi_bwt.allele_mapping_data.txt`: Read mapping at allele level
2. `{sample}_rgi_bwt.gene_mapping_data.txt`: Read mapping at gene level
3. `{sample}_rgi_bwt.sorted.length_100.bam`: BAM file with mapped reads
4. `{sample}_rgi_bwt.overall_mapping_stats.txt`: Overall statistics
5. `{sample}_rgi_bwt.reference_mapping_stats.txt`: Reference match statistics

### RGI K-mer Module (Pathogen-of-Origin)

**Module**: `modules/local/rgi_kmer/main.nf`

**Purpose**: Predict pathogen-of-origin for AMR genes using k-mer classification

**Inputs**:
- `tuple val(meta), path(bam)` - BAM file from RGI bwt
- `path(card_db)` - CARD database directory

**Process Configuration**:
- **Label**: `process_medium` (8 CPUs, 36GB RAM)
- **Container**: `quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0`

**Command Structure**:
```bash
rgi kmer_query \
    --bwt \
    --rgi_bwt_bam ${bam} \
    --kmer_size ${params.rgi_kmer_size} \
    --threads ${task.cpus} \
    --minimum 10 \
    --output ${prefix}_rgi_kmer \
    --local
```

**Outputs**:
- `tuple val(meta), path("*_61mer_analysis.json")` - K-mer analysis results
- `tuple val(meta), path("*_61mer_analysis.txt")` - Tab-delimited results
- `path("*.json")` - JSON output for downstream processing

**Output Files**:
1. `{sample}_rgi_kmer_61mer_analysis.json`: JSON format k-mer results
2. `{sample}_rgi_kmer_61mer_analysis.txt`: Tab-delimited k-mer results
3. `{sample}_rgi_kmer.allele.txt`: Allele-level pathogen predictions
4. `{sample}_rgi_kmer.gene.txt`: Gene-level pathogen predictions

### Report Generation Module

**Module**: `modules/local/rgi_report/main.nf`

**Purpose**: Aggregate RGI results across samples and generate summary reports

**Inputs**:
- `path(allele_files)` - All allele mapping files (collected)
- `path(gene_files)` - All gene mapping files (collected)
- `path(kmer_files)` - All k-mer analysis files (collected)

**Process Configuration**:
- **Label**: `process_low` (4 CPUs, 12GB RAM)
- **Container**: Custom Python container with pandas/matplotlib

**Outputs**:
- `path("RGI_summary_report.csv")` - Combined summary table
- `path("RGI_*.png")` - Visualization plots

**Report Contents**:
1. AMR gene family distribution across samples
2. Drug class resistance profiles
3. Resistance mechanism breakdown
4. Pathogen-of-origin predictions
5. Coverage and depth statistics

## Integration into Main Workflow

### Location in Pipeline

Insert RGI analysis after QC subworkflow, parallel to existing read-level ARG prediction:

```groovy
// In main.nf, after QC subworkflow (line ~241)

//
// RGI AMR PREDICTION IN READS
//
if ( params.rgi_prediction ) {
    // RGI bwt: Align reads to CARD
    ch_rgi_bwt = RGI_BWT(
        ch_clean_reads.combine(PREPARE_DATABASES.out.rgi_card_db)
    )
    
    // RGI kmer: Pathogen-of-origin prediction
    ch_rgi_kmer = RGI_KMER(
        ch_rgi_bwt.bam.combine(PREPARE_DATABASES.out.rgi_card_db)
    )
    
    // Generate summary report
    RGI_REPORT(
        ch_rgi_bwt.allele_mapping.collect(),
        ch_rgi_bwt.gene_mapping.collect(),
        ch_rgi_kmer.kmer_results.collect()
    )
}
```

### Channel Flow

```
INPUT_CHECK.out.reads
    ↓
QC.out.reads (ch_clean_reads)
    ↓
    ├─→ RGI_BWT(reads, card_db)
    │       ↓
    │       ├─→ allele_mapping_data.txt
    │       ├─→ gene_mapping_data.txt
    │       └─→ sorted.bam
    │               ↓
    │           RGI_KMER(bam, card_db)
    │               ↓
    │               └─→ kmer_analysis.json
    │
    └─→ (continue to existing TAXONOMY, ASSEMBLY, etc.)
```

## Configuration Updates

### nextflow.config

Add RGI parameters section:

```groovy
params {
    // RGI AMR prediction options
    rgi_prediction         = false         // Enable RGI analysis
    rgi_include_wildcard   = true          // Include WildCARD variants
    rgi_aligner            = 'kma'         // Aligner: kma, bowtie2, bwa
    rgi_kmer_size          = 61            // K-mer size for pathogen prediction
    rgi_min_kmer_coverage  = 10            // Minimum k-mer coverage
    
    // Custom RGI database paths
    custom_rgi_card_json   = null
    custom_rgi_wildcard    = null
}
```

### config/databases.config

Add RGI database preparation:

```groovy
process {
    withName: 'RGI_LOAD' {
        label = 'process_download'
        publishDir = [
            path: { "${params.databases_dir}/rgi" },
            mode: params.publish_dir_mode,
            pattern: '*'
        ]
    }
}
```

### config/modules.config

Add RGI module configurations:

```groovy
process {
    withName: 'RGI_BWT' {
        ext.args = [
            "--aligner ${params.rgi_aligner}",
            params.rgi_include_wildcard ? "--include_wildcard" : ""
        ].join(' ').trim()
        publishDir = [
            path: { "${params.output}/05_arg_prediction/read_level/rgi/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: '*.txt'
        ]
    }
    
    withName: 'RGI_KMER' {
        ext.args = [
            "--kmer_size ${params.rgi_kmer_size}",
            "--minimum ${params.rgi_min_kmer_coverage}"
        ].join(' ').trim()
        publishDir = [
            path: { "${params.output}/05_arg_prediction/read_level/rgi_kmer/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: '*.{json,txt}'
        ]
    }
    
    withName: 'RGI_REPORT' {
        publishDir = [
            [
                path: { "${params.output}/05_arg_prediction/read_level/rgi_summary" },
                mode: params.publish_dir_mode,
                pattern: '*.csv'
            ],
            [
                path: { "${params.output}/05_arg_prediction/read_level/rgi_summary" },
                mode: params.publish_dir_mode,
                pattern: '*.png'
            ]
        ]
    }
}
```

## Database Preparation Workflow

### PREPARE_DATABASES Subworkflow Update

Add to `subworkflows/local/prepare_databases.nf`:

```groovy
include { RGI_LOAD } from '../../modules/local/rgi_load/main'

workflow PREPARE_DATABASES {
    // ... existing code ...
    
    // RGI database preparation
    ch_rgi_card_db = Channel.empty()
    if ( params.rgi_prediction ) {
        if ( params.custom_rgi_card_json ) {
            // Use custom database
            ch_rgi_card_db = Channel.fromPath(params.custom_rgi_card_json)
        } else {
            // Download and prepare CARD database
            RGI_LOAD()
            ch_rgi_card_db = RGI_LOAD.out.card_db
        }
    }
    
    emit:
    // ... existing outputs ...
    rgi_card_db = ch_rgi_card_db
}
```

## Output Structure

```
results/
└── 05_arg_prediction/
    └── read_level/
        ├── rgi/
        │   ├── sample1/
        │   │   ├── sample1_rgi_bwt.allele_mapping_data.txt
        │   │   ├── sample1_rgi_bwt.gene_mapping_data.txt
        │   │   ├── sample1_rgi_bwt.overall_mapping_stats.txt
        │   │   └── sample1_rgi_bwt.reference_mapping_stats.txt
        │   └── sample2/
        │       └── ...
        ├── rgi_kmer/
        │   ├── sample1/
        │   │   ├── sample1_rgi_kmer_61mer_analysis.json
        │   │   ├── sample1_rgi_kmer_61mer_analysis.txt
        │   │   ├── sample1_rgi_kmer.allele.txt
        │   │   └── sample1_rgi_kmer.gene.txt
        │   └── sample2/
        │       └── ...
        └── rgi_summary/
            ├── RGI_summary_report.csv
            ├── RGI_amr_gene_family_distribution.png
            ├── RGI_drug_class_profile.png
            └── RGI_pathogen_origin_heatmap.png
```

## Implementation Steps

### Phase 1: Database Module
1. Create `modules/local/rgi_load/main.nf`
2. Implement CARD database download and preparation
3. Add WildCARD variant support
4. Test database preparation independently

### Phase 2: RGI BWT Module
1. Create `modules/local/rgi_bwt/main.nf`
2. Implement read alignment with KMA
3. Handle paired-end and singleton reads
4. Test with sample data

### Phase 3: RGI K-mer Module
1. Create `modules/local/rgi_kmer/main.nf`
2. Implement pathogen-of-origin prediction
3. Parse and format k-mer results
4. Test with RGI bwt output

### Phase 4: Report Module
1. Create `modules/local/rgi_report/main.nf`
2. Develop Python/R script for aggregation
3. Generate summary tables and plots
4. Test with multiple samples

### Phase 5: Integration
1. Update `subworkflows/local/prepare_databases.nf`
2. Update `main.nf` workflow
3. Add configuration parameters
4. Update help documentation

### Phase 6: Testing & Validation
1. Test with pipeline test data
2. Validate against known AMR samples
3. Compare with existing ARG prediction tools
4. Performance benchmarking

## Best Practices & Considerations

### Input Data Quality
- RGI bwt assumes pre-processed reads (QC, trimming, deduplication)
- BugBuster's FASTP QC step provides suitable input
- Consider read deduplication if not already implemented

### Database Selection
- **CARD only**: High specificity, clinical focus
- **CARD + WildCARD**: Better for environmental/non-clinical samples
- WildCARD increases allelic diversity but inflates allele network problem
- Recommend summarizing at AMR Gene Family level when using WildCARD

### Aligner Selection
- **KMA** (recommended): Best for redundant databases, resolves similar sequences
- **Bowtie2**: Faster but may align reads across multiple similar alleles
- **BWA**: Alternative, similar performance to Bowtie2

### Resource Requirements
- RGI bwt: ~8 CPUs, 36GB RAM per sample
- RGI kmer: ~8 CPUs, 36GB RAM per sample
- Database preparation: ~2 CPUs, 4GB RAM, ~20GB disk space
- WildCARD adds ~50GB disk space

### Performance Optimization
- Parallel processing of samples (Nextflow handles automatically)
- Use local database caching to avoid re-downloads
- Consider read downsampling for very large datasets

## Comparison with Existing Tools

### vs. KARGA/KARGVA
- **RGI**: Alignment-based, comprehensive CARD database, pathogen-of-origin
- **KARGA**: K-mer based, MEGARes database, faster but less comprehensive
- **Complementary**: Different databases and approaches provide validation

### vs. DeepARG
- **RGI**: Read-level analysis, no assembly required
- **DeepARG**: Protein-level, requires assembly and ORF prediction
- **Complementary**: Read-level (RGI) + contig-level (DeepARG) = comprehensive

### vs. ARG-OAP
- **RGI**: CARD database, detailed annotations, pathogen prediction
- **ARG-OAP**: SARG database, structured ARG types
- **Complementary**: Different databases provide cross-validation

## Documentation Updates

### README.md
Add RGI section to pipeline features and usage examples

### Help Text
Update `main.nf` help message with RGI parameters

### Example Command
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --quality_control true \
    --rgi_prediction true \
    --rgi_include_wildcard true \
    -profile docker
```

## Future Enhancements

1. **Mutation Screening**: Add support for protein variant models when RGI implements SNP screening
2. **Visualization**: Interactive HTML reports with AMR gene networks
3. **Integration**: Link RGI results with taxonomy data for pathogen validation
4. **Comparative Analysis**: Cross-reference RGI, KARGA, and DeepARG results
5. **Prevalence Data**: Incorporate CARD prevalence statistics into reports

## References

- RGI GitHub: https://github.com/arpcard/rgi
- CARD Database: https://card.mcmaster.ca/
- RGI Documentation: https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst
- K-mer Classifier: https://github.com/arpcard/rgi/blob/master/docs/rgi_kmer.rst
- Citation: Alcock et al. 2023. CARD 2023. Nucleic Acids Research, 51, D690-D699

---

**Document Version**: 1.0  
**Date**: 2026-01-25  
**Author**: BugBuster Development Team
