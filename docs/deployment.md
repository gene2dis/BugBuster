# BugBuster Deployment Guide

This guide covers deploying BugBuster on different compute infrastructures.

## Table of Contents

- [Local Execution](#local-execution)
- [HPC Clusters (SLURM)](#hpc-clusters-slurm)
- [AWS Batch](#aws-batch)
- [Google Cloud](#google-cloud)
- [Azure Batch](#azure-batch)
- [Seqera Platform](#seqera-platform)

---

## Local Execution

### Docker (Recommended)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker
```

### Singularity

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile singularity
```

---

## HPC Clusters (SLURM)

### Basic SLURM Execution

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output /scratch/user/results \
    -profile slurm,singularity
```

### With Custom Settings

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output /scratch/user/results \
    -profile slurm,singularity \
    --slurm_queue 'general' \
    --slurm_account 'my_allocation' \
    --singularity_cache '/shared/containers'
```

### Creating an Institutional Profile

1. Copy the template:
   ```bash
   cp conf/institutional.config conf/my_institution.config
   ```

2. Customize settings in the file

3. Add to `nextflow.config` profiles:
   ```groovy
   my_institution {
       includeConfig 'conf/my_institution.config'
   }
   ```

4. Run with your profile:
   ```bash
   nextflow run main.nf -profile my_institution
   ```

---

## AWS Batch

### Prerequisites

1. **AWS CLI configured** with credentials
2. **S3 bucket** for work directory and results
3. **AWS Batch** compute environment and job queue configured

### Environment Variables

```bash
export AWS_ACCESS_KEY_ID="your_access_key"
export AWS_SECRET_ACCESS_KEY="your_secret_key"
export AWS_DEFAULT_REGION="us-east-1"
```

### Basic Execution

```bash
nextflow run main.nf \
    --input s3://bucket/samplesheet.csv \
    --output s3://bucket/results \
    -profile aws,docker \
    --aws_queue 'my-batch-queue' \
    --aws_region 'us-east-1' \
    -work-dir s3://bucket/work
```

### With Spot Instances

AWS Batch automatically handles spot instance interruptions with the retry strategy configured in `conf/aws.config`.

### Fusion Filesystem (Recommended for Performance)

Enable Fusion for improved S3 performance (requires Nextflow 23.10+):

```bash
nextflow run main.nf \
    --input s3://bucket/samplesheet.csv \
    --output s3://bucket/results \
    -profile aws,docker \
    -with-wave \
    -with-fusion
```

---

## Google Cloud

### Prerequisites

1. **Google Cloud SDK** authenticated
2. **GCS bucket** for work directory and results
3. **Google Cloud Batch** or **Life Sciences API** enabled

### Authentication

```bash
gcloud auth application-default login
# Or set service account
export GOOGLE_APPLICATION_CREDENTIALS="/path/to/service-account.json"
```

### Basic Execution

```bash
nextflow run main.nf \
    --input gs://bucket/samplesheet.csv \
    --output gs://bucket/results \
    -profile gcp,docker \
    --gcp_project 'my-project-id' \
    --gcp_region 'us-central1' \
    -work-dir gs://bucket/work
```

### With Preemptible/Spot Instances

```bash
nextflow run main.nf \
    --input gs://bucket/samplesheet.csv \
    --output gs://bucket/results \
    -profile gcp,docker \
    --gcp_project 'my-project-id' \
    --gcp_spot true \
    -work-dir gs://bucket/work
```

---

## Azure Batch

### Prerequisites

1. **Azure Batch account** created
2. **Azure Storage account** for work directory
3. **Container registry** (optional, for custom images)

### Environment Variables

```bash
export AZURE_BATCH_ACCOUNT_NAME="mybatchaccount"
export AZURE_BATCH_ACCOUNT_KEY="batch_key"
export AZURE_STORAGE_ACCOUNT_NAME="mystorageaccount"
export AZURE_STORAGE_ACCOUNT_KEY="storage_key"
```

### Basic Execution

```bash
nextflow run main.nf \
    --input az://container/samplesheet.csv \
    --output az://container/results \
    -profile azure,docker \
    --azure_region 'eastus' \
    -work-dir az://container/work
```

### With Low-Priority (Spot) Instances

```bash
nextflow run main.nf \
    --input az://container/samplesheet.csv \
    --output az://container/results \
    -profile azure,docker \
    --azure_spot true \
    -work-dir az://container/work
```

---

## Seqera Platform

[Seqera Platform](https://seqera.io/) (formerly Nextflow Tower) provides a web interface for launching and monitoring pipelines.

### Setup

1. Create account at [cloud.seqera.io](https://cloud.seqera.io)
2. Generate access token
3. Configure compute environment

### Launch via CLI

```bash
export TOWER_ACCESS_TOKEN="your_token"

nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker \
    -with-tower
```

### Launch via Web Interface

1. Add pipeline to Launchpad
2. Configure compute environment
3. Set parameters
4. Launch run

---

## Resource Recommendations

### Minimum Requirements

| Component | Memory | CPUs | Storage |
|-----------|--------|------|---------|
| QC + Taxonomy | 16 GB | 4 | 50 GB |
| Assembly | 64 GB | 16 | 100 GB |
| Binning | 128 GB | 16 | 200 GB |
| Full Pipeline | 256 GB | 32 | 500 GB |

### Pre-downloaded Databases

For better performance, pre-download databases to shared storage:

```bash
# Kraken2 Standard-8 (~8 GB)
# Sourmash GTDB (~3 GB)
# CheckM2 (~3 GB)
# GTDB-TK (~85 GB)
```

Then specify paths:

```bash
nextflow run main.nf \
    --custom_kraken_db /shared/db/kraken2/standard-8 \
    --custom_checkm2_db /shared/db/checkm2/uniref100.KO.1.dmnd \
    --custom_gtdbtk_db /shared/db/gtdbtk/release220 \
    ...
```

---

## Troubleshooting

### Common Issues

1. **Out of memory**: Increase `--max_memory` or use larger instance types
2. **Spot instance preemption**: Pipeline will automatically retry
3. **Container pull failures**: Check network access or use pre-cached images
4. **S3/GCS access denied**: Verify IAM roles and bucket permissions

### Debug Mode

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --output ./results \
    -profile docker \
    -resume \
    -with-trace \
    -with-report \
    -with-dag
```

### Getting Help

- [BugBuster GitHub Issues](https://github.com/gene2dis/BugBuster/issues)
- [Nextflow Documentation](https://nextflow.io/docs/latest/)
- [nf-core Community](https://nf-co.re/)
