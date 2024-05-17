## Introducción

**Bacterial Unraveling and Genomic Binning with Up-Scale Throughput, Efficient and Reproducible**

**BUGBUSTER** is a bioinformatics best-practice analysis pipeline for microbial metagenomic and gene resistance analysis pipeline.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.

## Pipeline summary

1. Read QC, clean, and filter reads. [`FastP`](https://github.com/OpenGene/fastp)
2. Remove all samples that not have at least 1.000.000 reads after quality filter
3. Remove host contamintant reads. [`Bowtie2`](https://github.com/BenLangmead/bowtie2)
4. If requested Antibiotic resistance prediction at read level using KARGA and KARGVA [`KARGA`](https://github.com/DataIntellSystLab/KARGA), [`KARGVA`](https://github.com/DataIntellSystLab/KARGVA)
5. Taxonomic profile [`Kraken2`](https://ccb.jhu.edu/software/kraken2/)
6. Abundance estimation [`Bracken`](https://github.com/jenniferlu717/Bracken)
7. Read traceback and taxonomic reports
8. Genome assembly [`Megahit`](https://github.com/voutcn/megahit)
9. Contig filter [`BBmap`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
10. Binning [`Metabat2`](https://bitbucket.org/berkeleylab/metabat/src/master/), [`SemiBin`](https://github.com/BigDataBiology/SemiBin) and [`Autometa`](https://autometa.readthedocs.io/en/latest/getting-started.html)
11. Binning refinement [`MetaWrap`](https://github.com/bxlab/metaWRAP)
12. Bin quality prediction [`CheckM2`](https://github.com/chklovski/CheckM2)
13. Bin taxonomic prediction [`GTDB-TK`](https://github.com/Ecogenomics/GTDBTk)
14. Bin reports
15. If requested functional anotation of Reads, Contigs or/and Bins [`MetaCerberus`](https://github.com/raw-lab/MetaCerberus)
16. If requested **(work in progress)** ARG clustering and horizontal gene transfer inference [`mmseqs2`](https://github.com/soedinglab/MMseqs2), [`galaxy-tool-lca`](https://github.com/naturalis/galaxy-tool-lca)
17. Assembly modes: "coassembly", "assembly"

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) for full pipeline reproducibility.

Currently the pipeline has been tested with Docker.

3. Download the pipeline, either cloning the repository or downloading the zip file

## Running the pipeline

### Sample sheet file

First you need to create a Samplesheet file, that contain the name of the samples and the location of the reads. This is an example of this file:

```
sample,r1,r2,s
SRR9040400,/home/ffuentes/Raw_data/SRR9040400_1.fastq.gz,/home/ffuentes/Raw_data/SRR9040400_2.fastq.gz,/home/ffuentes/Raw_data/SRR9040400.fastq.gz
```

The file **always** has to include the header `sample,r1,r2,s`, "s" column can be empty 

# Pipeline parameters

The file `nextflow.config` contains all the parameteres used by the pipeline, including path to database files. Currently the path work in our server (_Arrakis_), but if you are running elsewhere, these need to be updated. 

#### Starting the pipeline
```
Usage:

The typical command for running the pipeline is as follows:

nextflow run main.nf --input "path/to/samples_sheet" --output "path/to/output" -resume (recomended)

Mandatory arguments:
   \color{green}--input                        Input csv file with: samples names, path of all fastq files, and optionaly singletons.
                                  colnames required: "sample,r1,r2,s" if don't have singletons colname "s" can be empty

   --output                       Path to output dir

Optional arguments:
   --assembly_mode                Mode of assembly, avaible options: "coassembly", "assembly"
                                  (default: assembly)
                                  coassembly: all samples are processing in one only data set
                                  assembly: all samples are processing individualy

   --single_assembly              If coassembly type is choosen, single assembly will aditionaly generate indivdual assembly for all samples
                                  (default: false)

   --arg_clustering               ARG gene prediction and clustering for horizontal gene transfer inference
                                  (default: false)

   --read_arg_prediction          ARG and ARGV gene prediction at read level using KARGA and KARGVA
                                  (default: false)

   --read_level_metacerberus      Read level functional annotation with metacerberus
                                  (default: false)

   --contig_level_metacerberus    Contig level functional annotation with metacebeus
                                  (default: false)

   --help                         Print this usage statement.
```
Additionally, all options can be modified in nextflow.config file
## Credits

gene2dis/BUGBUSTER was originally written by the Microbial Data Science Lab, Center for Bioinformatics and Integrative Biology, Universidad Andres Bello. Its development was led by Francisco A. Fuentes

We thank the following people for their extensive assistance in the development of this pipeline:

- Francisco A. Fuentes
- Juan A. Ugalde
- Shrek

![](https://pbs.twimg.com/media/DwoSzyiWkAA2ywL.jpg)
