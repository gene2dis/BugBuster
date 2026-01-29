# RGI Implementation Summary

## Overview

The RGI (Resistance Gene Identifier) AMR prediction functionality has been successfully implemented in the BugBuster pipeline. This implementation enables read-level AMR gene detection with pathogen-of-origin prediction using the CARD database.

## Implementation Date

**Completed**: January 25, 2026

## Files Created

### Modules

1. **`modules/local/rgi_load/main.nf`**
   - Downloads and prepares CARD database
   - Optionally includes WildCARD variants
   - Builds KMA indices for alignment

2. **`modules/local/rgi_bwt/main.nf`**
   - Aligns metagenomic reads to CARD AMR alleles
   - Uses KMA aligner (default) for optimal performance
   - Outputs allele-level and gene-level mapping results

3. **`modules/local/rgi_kmer/main.nf`**
   - K-mer based pathogen-of-origin prediction
   - Analyzes BAM output from RGI bwt
   - Identifies likely source organisms for detected AMR genes

4. **`modules/local/rgi_report/main.nf`**
   - Aggregates results across all samples
   - Generates summary tables and visualizations
   - Creates AMR gene family, drug class, and mechanism plots

## Files Modified

### Main Workflow

1. **`main.nf`**
   - Added RGI module imports (lines 173-176)
   - Added RGI prediction workflow section (lines 248-270)
   - Updated help text to include `--rgi_prediction` option
   - Updated run configuration logging

### Subworkflows

2. **`subworkflows/local/prepare_databases.nf`**
   - Added RGI_LOAD module import
   - Added RGI database preparation logic (lines 163-177)
   - Added `rgi_card_db` to emit block

### Configuration

3. **`nextflow.config`**
   - Added `rgi_prediction` feature toggle (line 148)
   - Added `custom_rgi_card_db` parameter (line 181)
   - Added RGI-specific parameters section (lines 263-268):
     - `rgi_card_version`
     - `rgi_include_wildcard`
     - `rgi_aligner`
     - `rgi_kmer_size`
     - `rgi_min_kmer_coverage`

4. **`config/modules.config`**
   - Added RGI module process configurations (lines 715-756)
   - Configured publishDir for all RGI processes
   - Set appropriate process labels

## Usage

### Basic Usage

Enable RGI prediction in your pipeline run:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    -profile docker
```

### With WildCARD Variants

Include extended allelic diversity (recommended for environmental samples):

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --rgi_include_wildcard true \
    -profile docker
```

### Using Existing CARD Database

If you have a pre-prepared CARD database:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --custom_rgi_card_db /path/to/card_database \
    -profile docker
```

### Custom Parameters

Adjust RGI-specific parameters:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --rgi_aligner kma \
    --rgi_kmer_size 61 \
    --rgi_min_kmer_coverage 10 \
    -profile docker
```

## Output Structure

```
results/
└── 05_arg_prediction/
    └── read_level/
        ├── rgi/
        │   └── {sample_id}/
        │       ├── {sample}_rgi_bwt.allele_mapping_data.txt
        │       ├── {sample}_rgi_bwt.gene_mapping_data.txt
        │       ├── {sample}_rgi_bwt.sorted.length_100.bam
        │       ├── {sample}_rgi_bwt.overall_mapping_stats.txt
        │       └── {sample}_rgi_bwt.reference_mapping_stats.txt
        ├── rgi_kmer/
        │   └── {sample_id}/
        │       ├── {sample}_rgi_kmer_61mer_analysis.json
        │       ├── {sample}_rgi_kmer_61mer_analysis.txt
        │       ├── {sample}_rgi_kmer.allele.txt (optional)
        │       └── {sample}_rgi_kmer.gene.txt (optional)
        └── rgi_summary/
            ├── RGI_summary_report.csv
            ├── RGI_amr_gene_family_distribution.png
            ├── RGI_drug_class_profile.png
            ├── RGI_resistance_mechanisms.png
            └── RGI_sample_amr_counts.png
```

## Key Features

### Database Flexibility

- **Automatic download**: Pipeline downloads CARD database automatically
- **Pre-prepared database**: Use existing CARD database with `--custom_rgi_card_db`
- **WildCARD support**: Optional extended allelic variants for non-clinical samples
- **Version control**: Specify CARD version or use latest

### Analysis Capabilities

- **Read-level AMR detection**: No assembly required
- **Pathogen-of-origin prediction**: K-mer based source organism identification
- **Comprehensive annotations**: Drug class, resistance mechanism, AMR gene family
- **Multi-sample analysis**: Aggregated reports across all samples

### Performance

- **Optimized aligner**: KMA recommended for redundant databases
- **Parallel processing**: Per-sample parallelization via Nextflow
- **Resource efficient**: ~8 CPUs, 36GB RAM per sample
- **Scalable**: Handles multiple samples efficiently

## Integration with Existing Tools

RGI complements existing ARG prediction tools in BugBuster:

- **KARGA/KARGVA**: Different database (MEGARes vs CARD), k-mer vs alignment
- **DeepARG**: Contig-level vs read-level, provides cross-validation
- **ARG-OAP**: Different database (SARG), structured ARG types

Running multiple tools provides comprehensive AMR gene detection and validation.

## Parameters Reference

### Feature Toggle

- `--rgi_prediction` (default: false) - Enable RGI AMR prediction

### Database Parameters

- `--rgi_card_version` (default: 'latest') - CARD database version
- `--rgi_include_wildcard` (default: true) - Include WildCARD variants
- `--custom_rgi_card_db` (default: null) - Path to pre-prepared CARD database

### Analysis Parameters

- `--rgi_aligner` (default: 'kma') - Read aligner: kma, bowtie2, bwa
- `--rgi_kmer_size` (default: 61) - K-mer size for pathogen prediction
- `--rgi_min_kmer_coverage` (default: 10) - Minimum k-mer coverage threshold

## Resource Requirements

### Per Sample

- **CPUs**: 8 cores
- **Memory**: 36 GB RAM
- **Time**: ~30-60 minutes (depends on read count and database size)

### Database Storage

- **CARD only**: ~500 MB
- **CARD + WildCARD**: ~20-50 GB
- **Temporary space**: ~100 GB during preparation

## Container Information

All RGI modules use the official biocontainer:

```
quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0
```

This container includes:
- RGI 6.0.3
- KMA 1.4.14
- Bowtie2 2.5.4
- BWA 0.7.18
- All required dependencies

## Testing

### Test Run

Test the implementation with a small dataset:

```bash
nextflow run main.nf \
    --input test_samplesheet.csv \
    --output ./test_results \
    --rgi_prediction true \
    --rgi_include_wildcard false \
    -profile docker
```

### Validation

Verify outputs:

1. Check RGI bwt outputs exist for each sample
2. Verify k-mer analysis completed
3. Review summary report and plots
4. Compare with existing ARG prediction results

## Troubleshooting

### Common Issues

**Issue**: Database download fails
- **Solution**: Check internet connection, try manual download (see implementation plan)

**Issue**: Out of memory errors
- **Solution**: Increase memory allocation or reduce sample size

**Issue**: No k-mer results
- **Solution**: Check BAM file quality, verify sufficient AMR gene coverage

**Issue**: WildCARD preparation takes too long
- **Solution**: Use CARD only for initial testing, prepare WildCARD separately

## Future Enhancements

Potential improvements for future versions:

1. **Mutation screening**: Add support for protein variant models when RGI implements SNP screening
2. **Interactive reports**: HTML reports with AMR gene networks
3. **Taxonomy integration**: Link RGI results with taxonomy data for validation
4. **Comparative analysis**: Cross-reference RGI, KARGA, and DeepARG results
5. **Prevalence data**: Incorporate CARD prevalence statistics

## Documentation

For detailed implementation information, see:

- **`docs/RGI_IMPLEMENTATION_PLAN.md`** - Complete implementation plan with manual database preparation
- **RGI GitHub**: https://github.com/arpcard/rgi
- **CARD Database**: https://card.mcmaster.ca/

## Citation

When using RGI in publications, please cite:

Alcock et al. 2023. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51, D690-D699. PMID: 36263822

## Support

For issues or questions:

1. Check the implementation plan documentation
2. Review RGI documentation: https://github.com/arpcard/rgi/tree/master/docs
3. Log an issue on the BugBuster GitHub repository
4. Contact CARD support: card@mcmaster.ca

---

**Implementation Status**: ✅ Complete and Ready for Testing

**Next Steps**:
1. Test with sample data
2. Validate outputs
3. Compare with existing ARG prediction tools
4. Document any issues or improvements needed
