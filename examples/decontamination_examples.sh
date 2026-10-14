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
  --quality_control true \
  -profile docker

# ==============================================================================
# Example 2: Using Pre-built Combined Index
# ==============================================================================
# Fastest option if you have a pre-built index
# Skips index building step

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --custom_decontamination_index /data/databases/contaminants_index \
  -profile docker

# ==============================================================================
# Example 3: Custom Host Genome (Mouse)
# ==============================================================================
# Build combined index with phiX and mouse genome

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --custom_host_fasta /data/genomes/mouse_GRCm39.fasta \
  -profile docker

# ==============================================================================
# Example 4: Multiple Contaminant Genomes
# ==============================================================================
# Remove reads mapping to phiX, human, mouse, and E. coli: --custom_host_fasta takes
# a comma-separated list, all built into one index. The bwa-mem3 build requests
# ~24 GB RAM per Gbp of genome (~140 GB here), so raise --max_memory

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --custom_phiX_fasta /data/genomes/phix.fasta \
  --custom_host_fasta /data/genomes/human.fasta,/data/genomes/mouse.fasta,/data/genomes/ecoli.fasta \
  --max_memory 150.GB \
  -profile docker

# ==============================================================================
# Example 5: Low Disk Space Mode
# ==============================================================================
# Enable automatic work-dir cleanup and clean-reads publishing
# Frees the work dir after a successful run; a completed low_disk run cannot be
# resumed (cleanup = true), an interrupted one can

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  -profile docker,low_disk

# ==============================================================================
# Example 6: Bowtie2 Fallback with Custom Parameters
# ==============================================================================
# The --bowtie_* options apply only under --read_aligner bowtie2 (the fallback
# for CPUs without AVX2 or nodes with < ~75 GB RAM); the default bwa-mem3 uses
# --decontam_min_score. bowtie2 runs in --local mode, so --bowtie_score_min takes
# the local form G,<a>,<b> (default G,15,6)

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --read_aligner bowtie2 \
  --bowtie_k 5 \
  --bowtie_score_min 'G,20,8' \
  -profile docker

# ==============================================================================
# Example 7: Skip Quality Control (No Decontamination)
# ==============================================================================
# Useful for pre-cleaned data or testing

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control false \
  -profile docker

# ==============================================================================
# Example 8: Building Combined Index Manually
# ==============================================================================
# Pre-build index for reuse across multiple runs

# Step 1: Concatenate FASTA files
cat phix.fasta human.fasta > contaminants.fasta

# Step 2: Build the bwa-mem3 index (default aligner; ~75 GB RAM for a human host).
# For --read_aligner bowtie2, build a Bowtie2 index instead (and add
# --read_aligner bowtie2 in step 3):
#   mkdir -p contaminants_index
#   bowtie2-build --threads 16 contaminants.fasta contaminants_index/contaminants
mkdir -p contaminants_index
bwa-mem3 index -t 16 -p contaminants_index/contaminants contaminants.fasta

# Step 3: Use in pipeline
nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --custom_decontamination_index contaminants_index \
  -profile docker

# ==============================================================================
# Example 9: Slurm Cluster with Custom Resources
# ==============================================================================
# Run on HPC with custom memory and CPU allocation

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  -profile slurm,singularity \
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
  -profile docker \
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
# Publish the downloaded databases and built indexes to a shared location. Later
# runs reuse them only when passed back with the --custom_* options (e.g.
# --custom_decontamination_index /shared/databases/bwamem3_index/contaminants_index)

nextflow run main.nf \
  --input samplesheet.csv \
  --output results \
  --quality_control true \
  --databases_dir /shared/databases \
  -profile docker

# ==============================================================================
# Performance Optimization Tips
# ==============================================================================

# Tip 1: Pre-build and cache the combined index
# Build once, use many times across different sample sets

# Tip 2: Publish clean reads for reuse outside the pipeline
# (publishing only - '-resume' remains the task cache)
nextflow run main.nf --input samplesheet.csv --output results --store_clean_reads true -profile docker

# Tip 3: Use the low_disk profile for constrained systems
# (enables Nextflow's 'cleanup = true' work-dir cleanup; a completed run cannot be resumed)
nextflow run main.nf --input samplesheet.csv --output results -profile docker,low_disk

# ==============================================================================
# Monitoring and Debugging
# ==============================================================================

# Check decontamination statistics
cat results/01_quality_control/summary/Reads_report.csv

# View pipeline execution timeline
open results/pipeline_info/execution_timeline_*.html

# Check resource usage
cat results/pipeline_info/execution_trace_*.txt

# View workflow DAG
open results/pipeline_info/pipeline_dag_*.svg
