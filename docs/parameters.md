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
8. [Assembly Parameters](#assembly-parameters)
9. [Binning Parameters](#binning-parameters)
10. [ARG Prediction Parameters](#arg-prediction-parameters)
11. [Functional Annotation Parameters](#functional-annotation-parameters)
12. [Alignment Parameters](#alignment-parameters)
13. [Clustering Parameters](#clustering-parameters)
14. [Resource Limit Parameters](#resource-limit-parameters)
15. [Advanced Parameters](#advanced-parameters)

---

## Input/Output Parameters

### `--input`
- **Type**: String (file path)
- **Required**: Yes
- **Description**: Path to CSV samplesheet containing sample information
- **Format**: CSV file with columns: `sample`, `r1`, `r2`, `s` (optional)
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
- **Options**: `kraken2`, `sourmash`, `none`
- **Description**: Tool for taxonomic profiling at read level
- **Example**: `--taxonomic_profiler kraken2`

### `--include_binning`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable metagenomic binning and refinement steps
- **Example**: `--include_binning true`

### `--read_arg_prediction`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable ARG and ARGV gene prediction at read level using KARGA and KARGVA
- **Example**: `--read_arg_prediction true`

### `--contig_tax_and_arg`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable contig-level taxonomy and ARG prediction using BlobTools and DeepARG
- **Example**: `--contig_tax_and_arg true`

### `--contig_level_metacerberus`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Enable contig-level functional annotation with MetaCerberus
- **Example**: `--contig_level_metacerberus true`

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
- **Description**: PhiX genome index selection for automatic download
- **Size**: 8.1 MB
- **Example**: `--phiX_index phiX174`

### `--host_db`
- **Type**: String
- **Default**: `human`
- **Description**: Host genome database for read filtering
- **Size**: 4.1 GB (CHM13 plus Y)
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

### `--checkm2_db`
- **Type**: String
- **Default**: `v3`
- **Description**: CheckM2 database version
- **Size**: 2.9 GB
- **Example**: `--checkm2_db v3`

### `--gtdbtk_db`
- **Type**: String
- **Default**: `release_220`
- **Description**: GTDB-TK database release version
- **Size**: 109 GB
- **Example**: `--gtdbtk_db release_220`

### `--databases_dir`
- **Type**: String (directory path)
- **Default**: `<output>/../databases`
- **Description**: Directory for storing downloaded databases (separate from results output)
- **Example**: `--databases_dir /shared/databases/bugbuster`

---

## Custom Database Paths

Override automatic downloads by providing custom database paths:

### `--custom_phiX_index`
- **Type**: String (directory path)
- **Description**: Path to custom PhiX Bowtie2 index directory
- **Example**: `--custom_phiX_index /path/to/phiX_index`

### `--custom_bowtie_host_index`
- **Type**: String (directory path)
- **Description**: Path to custom host Bowtie2 index directory
- **Example**: `--custom_bowtie_host_index /path/to/host_index`

### `--custom_kraken_db`
- **Type**: String (directory path)
- **Description**: Path to custom Kraken2 database directory
- **Example**: `--custom_kraken_db /path/to/kraken2_db`

### `--custom_sourmash_db`
- **Type**: Array of strings
- **Description**: Paths to custom Sourmash database files [kmer_file, lineages_file]
- **Example**: `--custom_sourmash_db '["/path/to/kmers.zip", "/path/to/lineages.csv"]'`

### `--custom_checkm2_db`
- **Type**: String (file path)
- **Description**: Path to custom CheckM2 database file
- **Example**: `--custom_checkm2_db /path/to/checkm2.dmnd`

### `--custom_gtdbtk_db`
- **Type**: String (directory path)
- **Description**: Path to custom GTDB-TK database directory
- **Example**: `--custom_gtdbtk_db /path/to/gtdbtk_r220`

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

### `--kraken_db_used`
- **Type**: String
- **Default**: `gtdb_release207`
- **Description**: Kraken2 database name for reports (metadata only)
- **Example**: `--kraken_db_used gtdb_release220`

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

### `--sourmash_db_name`
- **Type**: String
- **Default**: `gtdb_release_220`
- **Description**: Sourmash database name for reports (metadata only)
- **Example**: `--sourmash_db_name gtdb_release_220`

### `--sourmash_tax_rank`
- **Type**: String
- **Default**: `species`
- **Options**: `genus`, `species`, `strain`
- **Description**: Taxonomic rank for Sourmash classification
- **Example**: `--sourmash_tax_rank genus`

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
- **Description**: Number of top taxa to display in plots
- **Example**: `--taxonomy_top_n_taxa 20`

### `--create_phyloseq_rds`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Generate R phyloseq RDS files for downstream analysis
- **Example**: `--create_phyloseq_rds true`

---

## Assembly Parameters

### `--bbmap_lenght`
- **Type**: Integer
- **Default**: `1000`
- **Minimum**: `0`
- **Description**: Minimum contig length after BBMap filtering
- **Example**: `--bbmap_lenght 1500`

---

## Binning Parameters

### Basic Binning Parameters

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

---

## Functional Annotation Parameters

### `--metacerberus_hmm`
- **Type**: String
- **Default**: `"KOFam_all, COG, VOG, PHROG, CAZy"`
- **Options**: `KOFam_all`, `KOFam_eukaryote`, `KOFam_prokaryote`, `COG`, `VOG`, `PHROG`, `CAZy`
- **Description**: Comma-separated list of HMM databases to use for MetaCerberus
- **Example**: `--metacerberus_hmm "KOFam_prokaryote, COG, CAZy"`

### `--metacerberus_minscore`
- **Type**: Integer
- **Default**: `25`
- **Description**: Minimum HMM score for MetaCerberus
- **Example**: `--metacerberus_minscore 30`

### `--metacerberus_evalue`
- **Type**: Number
- **Default**: `1e-09`
- **Description**: Maximum E-value for MetaCerberus
- **Example**: `--metacerberus_evalue 1e-10`

---

## Alignment Parameters

Bowtie2 parameters for read alignment during host filtering:

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

## Clustering Parameters

MMseqs2 parameters for ARG clustering (used when `--arg_bin_clustering=true`):

### `--mmseqs_start_sens`
- **Type**: Integer
- **Default**: `2`
- **Description**: Starting sensitivity for MMseqs2
- **Example**: `--mmseqs_start_sens 3`

### `--mmseqs_s`
- **Type**: Integer
- **Default**: `7`
- **Description**: Sensitivity level for MMseqs2
- **Example**: `--mmseqs_s 8`

### `--mmseqs_sens_steps`
- **Type**: Integer
- **Default**: `3`
- **Description**: Number of sensitivity steps for MMseqs2
- **Example**: `--mmseqs_sens_steps 4`

### `--mmseqs_min_seq_id`
- **Type**: Number
- **Default**: `0.8`
- **Range**: 0-1
- **Description**: Minimum sequence identity for MMseqs2
- **Example**: `--mmseqs_min_seq_id 0.9`

### `--mmseqs_c`
- **Type**: Number
- **Default**: `0.7`
- **Range**: 0-1
- **Description**: Coverage threshold for MMseqs2
- **Example**: `--mmseqs_c 0.8`

### `--mmseqs_cov_mode`
- **Type**: Integer
- **Default**: `2`
- **Description**: Coverage mode for MMseqs2
- **Example**: `--mmseqs_cov_mode 1`

### `--mmseqs_e`
- **Type**: String
- **Default**: `1e-20`
- **Description**: E-value threshold for MMseqs2
- **Example**: `--mmseqs_e 1e-25`

### `--mmseqs_format_mode`
- **Type**: Integer
- **Default**: `4`
- **Description**: Output format mode for MMseqs2
- **Example**: `--mmseqs_format_mode 3`

### `--mmseqs_alignment_mode`
- **Type**: Integer
- **Default**: `3`
- **Description**: Alignment mode for MMseqs2
- **Example**: `--mmseqs_alignment_mode 2`

### `--mmseqs_max_seqs`
- **Type**: Integer
- **Default**: `10000`
- **Description**: Maximum number of sequences for MMseqs2
- **Example**: `--mmseqs_max_seqs 15000`

### `--mmseqs_format_output`
- **Type**: String
- **Default**: `empty,query,target,evalue,pident,qcov,tcov,tseq`
- **Description**: Output format fields for MMseqs2
- **Example**: `--mmseqs_format_output "query,target,pident,evalue"`

---

## Resource Limit Parameters

### `--max_cpus`
- **Type**: Integer
- **Default**: `16`
- **Description**: Maximum number of CPUs that can be requested for any single job
- **Example**: `--max_cpus 32`

### `--max_memory`
- **Type**: String
- **Default**: `128.GB`
- **Description**: Maximum amount of memory that can be requested for any single job
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

### `--validationShowHiddenParams`
- **Type**: Boolean
- **Default**: `false`
- **Description**: Show all parameters when using `--help`
- **Example**: `--validationShowHiddenParams`

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
    --custom_gtdbtk_db /shared/db/gtdbtk_r220 \
    --databases_dir /shared/databases \
    -profile docker
```

---

*BugBuster v1.0.0 - Complete Parameter Reference*
