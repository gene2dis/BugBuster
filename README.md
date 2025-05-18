![image](https://github.com/user-attachments/assets/a10c01f6-ef6c-40c4-a4ac-26a0c4f87564)
## Introduction

**Bacterial Unraveling and metaGenomic Binning with Up-Scale Throughput, Efficient and Reproducible** 

**BugBuster** is a bioinformatics best-practice analysis pipeline for microbial metagenomic and gene resistance.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.

## Pipeline summary
![Diagrama_BugBuster drawio](https://github.com/user-attachments/assets/40d02c04-e84e-48b4-b517-c93dccd68abc)


1. Read QC, clean, and filter reads. [`FastP`](https://github.com/OpenGene/fastp)
2. Remove all samples that not have at least 10.000.000 reads after quality filter. **Can be modified by the user**
3. Remove host contamintant reads. [`Bowtie2`](https://github.com/BenLangmead/bowtie2)
4. If requested Antibiotic resistance prediction at read level using KARGA and KARGVA [`KARGA`](https://github.com/DataIntellSystLab/KARGA), [`KARGVA`](https://github.com/DataIntellSystLab/KARGVA)
5. Normalization of predicted genes by estimating cell number with ARGs-OAP. [`ARGs-OAP`](https://github.com/xinehc/args_oap)
6. Taxonomic profile [`Kraken2`](https://ccb.jhu.edu/software/kraken2/) or [`Sourmash`](https://sourmash.readthedocs.io/en/latest/index.html)
7. Abundance estimation [`Bracken`](https://github.com/jenniferlu717/Bracken)
8. Unification of the results with Kraken-Biom and change of format to a Phyloseq object. [`Kraken-Biom`](https://github.com/smdabdoub/kraken-biom)
9. Read traceback and taxonomic reports.
10. Genome assembly [`Megahit`](https://github.com/voutcn/megahit)
11. Contig filter [`BBmap`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
12. Taxonomic annotation of contigs using Blastn and BlobTools. [`BlobTools`](https://github.com/DRL/blobtools), [`Blast`](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
13. Functional assignation of contigs with MetaCerberus. [`MetaCerberus`](https://github.com/raw-lab/MetaCerberus)
14. ORF prediction in contigs with Prodigal. [`Prodigal`](https://github.com/hyattpd/Prodigal)
15. Prediction of resistance genes at the contig level with DeepARG. [`DeepARG`](https://github.com/gaarangoa/deeparg)
16. Contig reports, scatter plot of taxonomy at Phylum level and scatter plot of resistance genes in contigs.
17. Binning [`Metabat2`](https://bitbucket.org/berkeleylab/metabat/src/master/), [`SemiBin`](https://github.com/BigDataBiology/SemiBin) and [`COMEBin`](https://github.com/ziyewang/COMEBin)
18. Binning refinement [`MetaWrap`](https://github.com/bxlab/metaWRAP)
19. Bin quality prediction [`CheckM2`](https://github.com/chklovski/CheckM2)
20. Bin taxonomic prediction [`GTDB-TK`](https://github.com/Ecogenomics/GTDBTk)
21. Bin reports.
22. If requested functional anotation of Bins **(work in progress)** [`MetaCerberus`](https://github.com/raw-lab/MetaCerberus)
23. If requested ARG clustering [`mmseqs2`](https://github.com/soedinglab/MMseqs2)
24. Assembly modes: "coassembly", "assembly", "none"

## Quick Start

Currently the pipeline has been tested with Docker.

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) for full pipeline reproducibility.

For install docker in Ubuntu, follow the instructions:
```bash
# Run the following command to uninstall all conflicting packages:
for pkg in docker.io docker-doc docker-compose docker-compose-v2 podman-docker containerd runc; do sudo apt-get remove $pkg; done

# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "${UBUNTU_CODENAME:-$VERSION_CODENAME}") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

# Install the docker packages
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# Verify that the installation is successful by running the hello-world image:
sudo docker run hello-world
```

**Once docker is installed and the pipeline is used for the first time, it will automatically download all the containers required for the execution of the pipeline with the desired configuration.**

3. Download the pipeline, either cloning the repository or downloading the zip file

4. Download databases and modify their paths in 'nextflow.config' file or use one of the default download databases from config/databases.config and modify the variable in the nextflow.config file.

## Databases

All databases can be automatically download in the first use of the pipeline and their paths will be stored as symbolic links in output_path/downloaded_db folder. The descriptions of the automatic download databases are in the config/databases.config file. However, you can use your own databases by writing the absolute paths in variables prefixed with custom_ in the nextflow.config file. 

**Note: Only the required databases for the requested tasks will be automatically download**

**Bowtie2:** must be directories with the genomes indexed with bowtie2 format
1. phiX_index = Default download from [`phage phiX174 geonome`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=genbank) (8.1 MB).
2. host_db = Default download from [`CHM13 plus Y bowtie2 index`](https://benlangmead.github.io/aws-indexes/bowtie) (4.1 GB).

**Taxonomic profiling:**
    
3A. kraken2_db = User can choice between Standard-8 (7.5 GB) and GTDB release 220 (497 GB) for automatic download. From [`kraken2 index`](https://benlangmead.github.io/aws-indexes/k2).
3B. sourmash_db = GTDB release 220 (17 GB) it's the only automatic default download in this instance (for now 👷). From [`sourmash kmers`](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs220/). 

**BBlobTools:**

4. taxdump_files = Default download from [`taxdump.tar.gz`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) (448 MB).

**BlastDB for contig taxonomic annotation:**

5. blast_db = Default download from [`ncbi nt`](https://ftp.ncbi.nlm.nih.gov/blast/db/) (434 GB).
6. deeparg_db = Default download from [`deeparg`](https://github.com/gaarangoa/deeparg) (4.8 GB).

**KARGA:**

7. karga_db = Default download from [`megares`](https://www.meglab.org/megares/download/) (9.2 MB).

**KARGVA:**

8. kargva_db = Default download from [`kargva_db`](https://github.com/DataIntellSystLab/KARGVA/tree/main) (1.5 MB).

**CheckM2**

9. checkm2_db = Default download from [`Checkm2_docs`](https://github.com/chklovski/CheckM2) (2.9 GB).

**GTDB-TK**

10. gtdbtk_db = Default download release 220 from [`gtdbtk_db`](https://ecogenomics.github.io/GTDBTk/installing/index.html) (109 GB).

## Running the pipeline

### Sample sheet file

First you need to create a Samplesheet file, that contain the name of the samples and the location of the reads. This is an example of this file:

```
sample,r1,r2,s
SRR9040400,/home/ffuentes/Raw_data/SRR9040400_1.fastq.gz,/home/ffuentes/Raw_data/SRR9040400_2.fastq.gz,/home/ffuentes/Raw_data/SRR9040400.fastq.gz
```

The file **always** has to include the header `sample,r1,r2,s` but also "s" column can be empty 

# Pipeline parameters

The file `nextflow.config` contains all the parameteres used by the pipeline.

#### Starting the pipeline

```
Usage:

The typical command for running the pipeline is as follows:

nextflow run main.nf --input "path/to/samples_sheet" --output "path/to/output" -profile local_docker -resume (recomended)

Mandatory arguments:
   --input                        Input csv file with: samples names, path of all fastq files, and optionaly singletons.
                                  colnames required: "sample,r1,r2,s" if don't have singletons colname "s" can be empty

   --output                       Path to output dir

Optional arguments:
   --quality_control              Include the reads filtering steps.
                                  (default: true)

   --taxonomic_profiler           Software for taxonomic classification in reads, avaible options: "kraken2", "sourmash", "none
                                  (default: kraken2)

   --assembly_mode                Mode of assembly, avaible options: "coassembly", "assembly", "none"
                                  (default: assembly)
                                  coassembly: all samples are processing in one only data set
                                  assembly: all samples are processing individualy

   --read_arg_prediction          ARG and ARGV gene prediction at read level using KARGA and KARGVA
                                  (default: false)

   --contig_level_metacerberus    Contig level functional annotation with metacebeus
                                  (default: false)

   --contig_tax_and_arg           Contig taxonomy and ARG prediction
                                  (default: false)

   --include_binning              Include the binning and binning refining steps. (default: false)

   --arg_bin_clustering           ARG gene prediction and clustering for horizontal gene transfer inference (WIP)
                                  (default: false)

   --help                         Print this usage statement.
```
Additionally, all options can be modified in nextflow.config file
## Credits

gene2dis/BUGBUSTER was originally written by the Microbial Data Science Lab, Center for Bioinformatics and Integrative Biology, Universidad Andres Bello. Its development was led by Francisco A. Fuentes

We thank the following people for their extensive assistance in the development of this pipeline:

- Francisco A. Fuentes
- Juan A. Ugalde
