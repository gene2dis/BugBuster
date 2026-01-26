# Using Pre-Downloaded WildCARD Database with RGI

This guide explains how to use an already downloaded WildCARD database with the BugBuster pipeline's RGI implementation.

## Overview

The pipeline supports three scenarios for RGI database usage:

1. **Automatic download** - Pipeline downloads both CARD and WildCARD
2. **Pre-prepared complete database** - Use a database that already includes WildCARD
3. **Separate CARD and WildCARD databases** - Combine them during pipeline execution

## Scenario 1: Automatic Download (Default)

The pipeline automatically downloads and prepares both CARD and WildCARD:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --rgi_include_wildcard true \
    -profile docker
```

**What happens:**
- Downloads CARD JSON database
- Downloads WildCARD variants
- Processes WildCARD annotations (30-60 minutes)
- Combines both into a single database
- Builds KMA indices

## Scenario 2: Pre-Prepared Complete Database

If you have a database that already includes both CARD and WildCARD:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --custom_rgi_card_db /path/to/complete_card_database \
    -profile docker
```

**Requirements:**
- Database must be RGI-loaded (created with `rgi load`)
- Must contain KMA indices
- WildCARD must already be integrated

**Database structure:**
```
complete_card_database/
├── card.json
├── card_database_v*.fasta
├── card_database_v*.fasta.comp.b
├── card_database_v*.fasta.length.b
├── card_database_v*.fasta.name
├── card_database_v*.fasta.seq.b
├── card_database_v*.fasta.index.b
├── wildcard_database_v*.fasta
├── wildcard_database_v*.fasta.comp.b
└── ... (other WildCARD KMA indices)
```

## Scenario 3: Separate CARD and WildCARD Databases (NEW)

If you have CARD and WildCARD databases stored separately, the pipeline can combine them:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --custom_rgi_card_db /path/to/card_only_database \
    --custom_rgi_wildcard /path/to/wildcard_directory \
    -profile docker
```

**What happens:**
- Uses your existing CARD database
- Processes WildCARD annotations from your directory
- Combines them into a new database
- Builds necessary KMA indices

### CARD Database Requirements

Your CARD-only database must contain:
```
card_only_database/
├── card.json
├── card_database_v*.fasta
└── localDB/ (RGI local database files)
```

### WildCARD Directory Requirements

Your WildCARD directory must contain the raw downloaded files:
```
wildcard_directory/
├── index-for-model-sequences.txt
├── nucleotide_fasta_protein_homolog_model_variants.fasta
├── nucleotide_fasta_protein_overexpression_model_variants.fasta
├── nucleotide_fasta_protein_variant_model_variants.fasta
├── nucleotide_fasta_rRNA_gene_variant_model_variants.fasta
└── ... (other variant files)
```

**Note:** Files can be gzipped (.gz) - the pipeline will decompress them.

## How to Download WildCARD Separately

If you only have CARD and need to download WildCARD:

```bash
# Create WildCARD directory
mkdir -p /path/to/wildcard_directory
cd /path/to/wildcard_directory

# Download WildCARD variants
wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants

# Extract files
tar -xjf wildcard_data.tar.bz2
rm wildcard_data.tar.bz2

# Decompress if needed
gunzip *.gz
```

## Performance Considerations

### Processing Time

| Scenario | CARD Download | WildCARD Processing | Total Time |
|----------|---------------|---------------------|------------|
| Automatic (CARD only) | ~5 min | - | ~5 min |
| Automatic (CARD + WildCARD) | ~10 min | 30-60 min | 40-70 min |
| Pre-prepared complete | - | - | Instant |
| Separate databases | - | 30-60 min | 30-60 min |

### Storage Requirements

- **CARD only**: ~500 MB
- **WildCARD raw files**: ~5 GB
- **CARD + WildCARD combined**: ~20-50 GB (with indices)

## Complete Examples

### Example 1: First-Time Setup with Automatic Download

```bash
# Let pipeline download everything
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --rgi_include_wildcard true \
    --databases_dir /shared/databases \
    -profile docker
```

Result: Database saved to `/shared/databases/rgi/` for future use.

### Example 2: Reuse Complete Database

```bash
# Use previously prepared database
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --custom_rgi_card_db /shared/databases/rgi/card_db \
    -profile docker
```

### Example 3: Add WildCARD to Existing CARD Database

```bash
# You have CARD at: /data/card_database
# You have WildCARD at: /data/wildcard

nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --custom_rgi_card_db /data/card_database \
    --custom_rgi_wildcard /data/wildcard \
    -profile docker
```

This creates a new combined database in the work directory and uses it for analysis.

### Example 4: CARD Only (No WildCARD)

```bash
# Use CARD without WildCARD variants
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --rgi_include_wildcard false \
    -profile docker
```

Or with custom CARD-only database:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    --rgi_prediction true \
    --custom_rgi_card_db /data/card_only \
    -profile docker
```

## Troubleshooting

### Issue: WildCARD Processing Fails

**Symptom**: Error during `rgi wildcard_annotation`

**Solutions:**
1. Ensure WildCARD files are decompressed
2. Check CARD and WildCARD versions match
3. Verify `index-for-model-sequences.txt` exists

### Issue: Version Mismatch

**Symptom**: "Version mismatch between CARD and WildCARD"

**Solution:**
```bash
# Check CARD version
grep '"version"' /path/to/card_database/card.json

# Check WildCARD version
grep "version" /path/to/wildcard/index-for-model-sequences.txt

# Download matching versions
```

### Issue: Insufficient Disk Space

**Symptom**: Pipeline fails during WildCARD processing

**Solution:**
- Ensure at least 100 GB free space in work directory
- Use `--databases_dir` to specify location with sufficient space
- Consider using CARD-only mode if space is limited

## Best Practices

1. **For Production**: Pre-prepare databases once and reuse with `--custom_rgi_card_db`
2. **For Development**: Use automatic download with `--rgi_include_wildcard false` for faster testing
3. **For Shared Systems**: Prepare databases in shared location, use custom paths
4. **For Environmental Samples**: Always include WildCARD for better coverage
5. **For Clinical Samples**: CARD-only may be sufficient

## Related Documentation

- **Full implementation details**: [`RGI_IMPLEMENTATION_PLAN.md`](RGI_IMPLEMENTATION_PLAN.md)
- **Manual database preparation**: See RGI_IMPLEMENTATION_PLAN.md Steps 1-5
- **Parameter reference**: [`parameters.md`](parameters.md)
- **Pipeline manual**: [`manual.md`](manual.md)

## Summary

| Parameter | Purpose | When to Use |
|-----------|---------|-------------|
| `--rgi_prediction` | Enable RGI | Always required for RGI analysis |
| `--rgi_include_wildcard` | Download WildCARD | Automatic download mode only |
| `--custom_rgi_card_db` | Use existing CARD database | When you have pre-prepared CARD |
| `--custom_rgi_wildcard` | Add separate WildCARD | When CARD and WildCARD are separate |

**Key Point**: Use both `--custom_rgi_card_db` and `--custom_rgi_wildcard` together when you have separately downloaded databases that need to be combined.
