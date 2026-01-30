#!/bin/bash
# BugBuster Decontamination Examples
# Demonstrates various ways to use the optimized single-pass decontamination

# ==============================================================================
# Example 1: Basic Usage (Default Settings)
# ==============================================================================
# Uses automatic download of phiX174 and human CHM13 genomes
# Builds combined index automatically

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true

# ==============================================================================
# Example 2: Using Pre-built Combined Index
# ==============================================================================
# Fastest option if you have a pre-built index
# Skips index building step

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --custom_decontamination_index /data/databases/contaminants_index

# ==============================================================================
# Example 3: Custom Host Genome (Mouse)
# ==============================================================================
# Build combined index with phiX and mouse genome

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --host_db mouse \
  --custom_host_fasta /data/genomes/mouse_GRCm39.fasta

# ==============================================================================
# Example 4: Multiple Contaminant Genomes
# ==============================================================================
# Remove reads mapping to phiX, human, mouse, and E. coli

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --custom_phiX_fasta /data/genomes/phix.fasta \
  --custom_host_fasta /data/genomes/human.fasta,/data/genomes/mouse.fasta,/data/genomes/ecoli.fasta

# ==============================================================================
# Example 5: Low Disk Space Mode
# ==============================================================================
# Enable progressive cleanup and use storeDir for clean reads
# Minimizes disk usage during pipeline execution

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  -profile low_disk

# ==============================================================================
# Example 6: Custom Bowtie2 Parameters
# ==============================================================================
# Adjust alignment sensitivity for decontamination

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --bowtie_k 5 \
  --bowtie_score_min 'L,-0.6,-0.6'

# ==============================================================================
# Example 7: Skip Quality Control (No Decontamination)
# ==============================================================================
# Useful for pre-cleaned data or testing

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control false

# ==============================================================================
# Example 8: Building Combined Index Manually
# ==============================================================================
# Pre-build index for reuse across multiple runs

# Step 1: Concatenate FASTA files
cat phix.fasta human.fasta > contaminants.fasta

# Step 2: Build Bowtie2 index
mkdir -p contaminants_index
bowtie2-build --threads 16 contaminants.fasta contaminants_index/contaminants

# Step 3: Use in pipeline
nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --custom_decontamination_index contaminants_index

# ==============================================================================
# Example 9: Slurm Cluster with Custom Resources
# ==============================================================================
# Run on HPC with custom memory and CPU allocation

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  -profile slurm \
  --max_memory 256.GB \
  --max_cpus 32

# ==============================================================================
# Example 10: Resume Failed Run
# ==============================================================================
# Resume from last successful checkpoint
# Combined index will be reused if already built

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  -resume

# ==============================================================================
# Example 11: Docker Container Execution
# ==============================================================================
# Run with Docker for reproducibility

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  -profile docker

# ==============================================================================
# Example 12: Custom Database Directory
# ==============================================================================
# Store databases in a specific location for sharing across runs

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --databases_dir /shared/databases

# ==============================================================================
# Performance Optimization Tips
# ==============================================================================

# Tip 1: Pre-build and cache the combined index
# Build once, use many times across different sample sets

# Tip 2: Use storeDir for clean reads in production
# Enables immediate work directory cleanup
nextflow run main.nf --store_clean_reads true

# Tip 3: Enable work directory cleanup for disk space
# Automatically removes intermediate files
nextflow run main.nf --enable_work_cleanup true

# Tip 4: Use low_disk profile for constrained systems
nextflow run main.nf -profile low_disk

# ==============================================================================
# Monitoring and Debugging
# ==============================================================================

# Check decontamination statistics
cat results/reads_report/reads_summary_report.tsv

# View pipeline execution timeline
open results/pipeline_info/execution_timeline_*.html

# Check resource usage
cat results/pipeline_info/execution_trace_*.txt

# View workflow DAG
open results/pipeline_info/pipeline_dag_*.svg
