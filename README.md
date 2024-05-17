## Introducción

**Bacterial Unraveling and Genomic Binning with Up-Scale Throughput, Efficient and Reproducible**

**BUGBUSTER** is a bioinformatics best-practice analysis pipeline for microbial metagenomic analysis pipeline.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

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

![](https://pbs.twimg.com/media/DwoSzyiWkAA2ywL.jpg)
