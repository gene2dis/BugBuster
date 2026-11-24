# BugBuster Deployment Guide

This guide covers deploying BugBuster on different compute infrastructures.
The cloud/HPC profile parameters used below (`aws_*`, `gcp_*`, `azure_*`,
`slurm_*`) are formally documented in
[`parameters.md` § Cloud & HPC Profile Parameters](parameters.md#cloud--hpc-profile-parameters);
the chosen cloud profile needs a work directory in its cloud storage, given with its
`--*_workdir` parameter or with `-work-dir` (see that section).

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

With Apptainer, use `-profile slurm,apptainer`, with `slurm` first: `conf/slurm.config` turns the
singularity engine on, and only an `apptainer` profile applied after it turns that off again
(`apptainer,slurm` enables both engines and Nextflow refuses to start). `--singularity_cache` and
`--singularity_run_options` apply to Apptainer too. Singularity needs version 3.6 or later (or any
Apptainer): the pipeline passes `--env` to some tasks.

### Creating an Institutional Profile

`conf/institutional.config` is a **template only** — no `institutional` profile is defined in `nextflow.config`, so it must be either registered as a profile (steps below) or passed directly with `-c conf/my_institution.config`.

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
   nextflow run main.nf \
       --input samplesheet.csv \
       --output ./results \
       -profile my_institution
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
    --aws_workdir s3://bucket/work
```

### With Spot Instances

A task killed by a spot reclaim is retried under the pipeline-wide rule in `nextflow.config` (exit statuses 130–145 and 104, up to 2 retries). The `errorStrategy` and `maxRetries = 3` in `conf/aws.config` are currently overridden by that rule ([#31](https://github.com/gene2dis/BugBuster/issues/31)). Nextflow's own spot handling is set with `aws.batch.maxSpotAttempts` in a `-c` config.

### Fusion Filesystem (Recommended for Performance)

Enable Fusion for improved S3 performance:

```bash
nextflow run main.nf \
    --input s3://bucket/samplesheet.csv \
    --output s3://bucket/results \
    -profile aws,docker \
    --aws_workdir s3://bucket/work \
    -with-wave \
    -with-fusion
```

---

## Google Cloud

### Prerequisites

1. **Google Cloud SDK** authenticated
2. **GCS bucket** for work directory and results
3. **Google Cloud Batch** API enabled

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
    --gcp_workdir gs://bucket/work
```

### With Preemptible/Spot Instances

```bash
nextflow run main.nf \
    --input gs://bucket/samplesheet.csv \
    --output gs://bucket/results \
    -profile gcp,docker \
    --gcp_project 'my-project-id' \
    --gcp_spot true \
    --gcp_workdir gs://bucket/work
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
    --azure_workdir az://container/work
```

### With Low-Priority (Spot) Instances

```bash
nextflow run main.nf \
    --input az://container/samplesheet.csv \
    --output az://container/results \
    -profile azure,docker \
    --azure_spot true \
    --azure_workdir az://container/work
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

### Machine Size and Databases

What each set of options asks of a node (memory, CPUs, database disk) is in the
manual's [hardware table](manual.md#hardware-requirements). On a cluster or in the
cloud, set `--max_memory` and `--max_cpus` to the largest node or instance type the
jobs can get: every request is capped by them. Decontamination and MEGAHIT request
72 GB each on the first attempt. The human decontamination index (~75 GB RAM with
bwa-mem3) is built once per `--databases_dir`, so put `--databases_dir` on storage
every node can read: the first run downloads and builds the databases there and
later runs reuse them. On small nodes or CPUs without AVX2, use
`--read_aligner bowtie2` (~6 GB, about an hour). To download databases by hand
instead (for compute nodes without internet access), see the manual,
[section 8](manual.md#manual-database-download).

> **Note — functional annotation branches on cloud executors:**
> `--contig_level_functional` runs under the docker-based `aws`, `gcp` and
> `azure` profiles using the pipeline-built eggNOG-mapper image
> (`ghcr.io/gene2dis/bugbuster-eggnog-mapper`; see the container note in
> [`manual.md`](manual.md)); its eggNOG-mapper search peaked at ~80 GB RAM on
> a real 5-sample co-assembly, so size the instance type / queue for it. The
> MAG branch (`--mag_level_functional`, Bakta) runs on any engine, and so
> is the read-level branch (`--read_level_functional woltka`) — but note its
> Bowtie2 alignment against WoLr2 needs ≥ 68 GB RAM per task, so size the
> cloud instance type / queue accordingly. The SUPER-FOCUS read backend
> (`--read_level_functional superfocus`) also runs on any engine; with
> `--superfocus_aligner mmseqs2` each task builds a ~13 GB k-mer index of the
> database in memory. The HUMAnN backend (`--read_level_functional humann`)
> runs on any engine too; its MetaPhlAn prescreen loads a ~20 GB Bowtie2
> index, so budget tens of GB of RAM per task (`process_high`).

---

## Troubleshooting

Cluster and cloud problems (SLURM jobs, access denied, quotas, spot instances,
authentication) are in [`troubleshooting.md`](troubleshooting.md#cloud-execution-issues),
with the rest of the troubleshooting. Questions:
[BugBuster GitHub Issues](https://github.com/gene2dis/BugBuster/issues),
[Nextflow documentation](https://nextflow.io/docs/latest/).
