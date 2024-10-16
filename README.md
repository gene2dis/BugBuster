![image](https://github.com/user-attachments/assets/a10c01f6-ef6c-40c4-a4ac-26a0c4f87564)
## Introduction

**Bacterial Unraveling and metaGenomic Binning with Up-Scale Throughput, Efficient and Reproducible** 

**BugBuster** is a bioinformatics best-practice analysis pipeline for microbial metagenomic and gene resistance.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.

## Pipeline summary
![Diagrama_BugBuster drawio](https://github.com/user-attachments/assets/9e54f1e8-9f6b-4181-a075-8bd2f6ebe078)

1. Read QC, clean, and filter reads. [`FastP`](https://github.com/OpenGene/fastp)
2. Remove all samples that not have at least 10.000.000 reads after quality filter. **Can be modified by the user**
3. Remove host contamintant reads. [`Bowtie2`](https://github.com/BenLangmead/bowtie2)
4. If requested Antibiotic resistance prediction at read level using KARGA and KARGVA [`KARGA`](https://github.com/DataIntellSystLab/KARGA), [`KARGVA`](https://github.com/DataIntellSystLab/KARGVA)
5. Normalization of predicted genes by estimating cell number with ARGs-OAP. [`ARGs-OAP`](https://github.com/xinehc/args_oap)
6. Taxonomic profile [`Kraken2`](https://ccb.jhu.edu/software/kraken2/)
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

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) for full pipeline reproducibility.

Currently the pipeline has been tested with Docker.

3. Download the pipeline, either cloning the repository or downloading the zip file

4. Download databases and modify their paths in 'nextflow.config' file.

## Databases

All database paths must be modified in the configuration file.

**Bowtie2:** must be directories with the genomes indexed with bowtie2 format
1. human_db = recommended download: [`Chm13plusY`](https://genome-idx.s3.amazonaws.com/bt/chm13.draft_v1.0_plusY.zip)
2. phyX_db = recommended download: [`phiX`](https://www.ncbi.nlm.nih.gov/nuccore/J02482.1?report=fasta) **Fasta file must be indexed with bowtie2**

**Kraken2:**

3. k2_gtdb_db = recommended download: [`gtdb_release_207_for_kraken2`](http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release207/)

**BBlobTools:**

4. ncbi_nodes_dmp = recommended download: nodes.dmp from zipped dir taxdump.tar.gz: [`taxdump.tar.gz`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
5. ncbi_names_dmp = recommended download: names.dmp from the same zipped dir of nodes.dmp

**CheckM for Metawrap:**

6. metawrap_db = recommended download: [`checkM_db_dir`](https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz)

**BlastDB for contig taxonomic annotation:**

7. nt_db = recommended download: [`ncbi_nt`](https://ftp.ncbi.nlm.nih.gov/blast/db/)
8. deeparg_db = Use the downloand command from deeparg software, for more info visit [`deeparg github`](https://github.com/gaarangoa/deeparg)

**KARGA:**

9. karga_db = must be a fasta with ARG genes. recommended download: [`megares_db`](https://www.meglab.org/downloads/megares_v3.00/megares_database_v3.00.fasta)

**¡KARGVA database for ARGV genes it's included in KARGVA container!**

**AUTOMETA, This is obsolete until future updates.**

10. ncbi_db = follow the instructions from [`Autometa_db_docs`](https://autometa.readthedocs.io/en/latest/databases.html) and generate a single directory with all files.

**CheckM2**

11. checkm_db = follow the instructions from [`Checkm2_docs`](https://github.com/chklovski/CheckM2)

**GTDB-TK**

12. gtdbtk-db = recommended download release 220 from [`gtdbtk_db`](https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz)

## Running the pipeline

### Sample sheet file

First you need to create a Samplesheet file, that contain the name of the samples and the location of the reads. This is an example of this file:

```
sample,r1,r2,s
SRR9040400,/home/ffuentes/Raw_data/SRR9040400_1.fastq.gz,/home/ffuentes/Raw_data/SRR9040400_2.fastq.gz,/home/ffuentes/Raw_data/SRR9040400.fastq.gz
```

The file **always** has to include the header `sample,r1,r2,s` but also "s" column can be empty 

# Pipeline parameters

The file `nextflow.config` contains all the parameteres used by the pipeline, including path to database files. Currently the path work in our server (_Arrakis_), but if you are running elsewhere, these need to be updated. 

#### Starting the pipeline

```
Usage:

The typical command for running the pipeline is as follows:

nextflow run main.nf --input "path/to/samples_sheet" --output "path/to/output" -resume (recomended)

Mandatory arguments:
   --input                        Input csv file with: samples names, path of all fastq files, and optionaly singletons.
                                  colnames required: "sample,r1,r2,s" if don't have singletons colname "s" can be empty

   --output                       Path to output dir

Optional arguments:
   --assembly_mode                Mode of assembly, avaible options: "coassembly", "assembly", "none"
                                  (default: assembly)
                                  coassembly: all samples are processing in one only data set
                                  assembly: all samples are processing individualy

   --single_assembly              If coassembly type is choosen, single assembly will aditionaly generate indivdual assembly for all samples
                                  (default: false)

   --read_arg_prediction          ARG and ARGV gene prediction at read level using KARGA and KARGVA
                                  (default: false)

   --contig_level_metacerberus    Contig level functional annotation with metacebeus
                                  (default: false)

   --contig_tax_and_arg           Contig taxonomy and ARG prediction
                                  (default: false)

   --arg_clustering               ARG gene prediction and clustering for horizontal gene transfer inference (WIP)
                                  (default: false)

   --help                         Print this usage statement.
```
Additionally, all options can be modified in nextflow.config file
## Credits

gene2dis/BUGBUSTER was originally written by the Microbial Data Science Lab, Center for Bioinformatics and Integrative Biology, Universidad Andres Bello. Its development was led by Francisco A. Fuentes

We thank the following people for their extensive assistance in the development of this pipeline:

- Francisco A. Fuentes
- Juan A. Ugalde
