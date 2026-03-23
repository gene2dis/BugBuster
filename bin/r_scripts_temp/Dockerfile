FROM rocker/r-ver:4.3.2
LABEL maintainer="ffuentessantander@gmail.com" 
RUN apt-get update && apt-get install zlib1g-dev && apt-get -y install libglpk-dev
RUN R -e "install.packages('igraph')"
RUN R -e "install.packages('tidyverse')"
RUN R -e "install.packages('plyr')"
RUN R -e "install.packages('BiocManager')" 
RUN R -e "BiocManager::install(c('phyloseq', 'microbiome', 'ComplexHeatmap'), update = FALSE)"
RUN R -e "install.packages('microViz',repos = c(davidbarnett = 'https://david-barnett.r-universe.dev', getOption('repos')))"
RUN R -e "install.packages('seqinr')"
RUN R -e "install.packages('rcartocolor')"
RUN R -e "install.packages('cowplot')"
RUN R -e "install.packages('optparse')"
RUN R -e "BiocManager::install('Rsubread')"
RUN R -e "install.packages('ape')"
RUN R -e "install.packages('RColorBrewer')"
RUN R -e "install.packages('reshape2')"
RUN R -e "BiocManager::install('edgeR')"
RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "install.packages('RSQLite')"
RUN R -e "install.packages('plyr')"
RUN R -e "BiocManager::install('topGO')"
RUN R -e "BiocManager::install('GO.db')"
RUN R -e "BiocManager::install('EnhancedVolcano')"
COPY Report_unify.R /mnt/
COPY Tax_unify_report.R /mnt/
COPY Bin_checkm_per_sample_plot.R /mnt/
COPY Bins_tax.R /mnt/
COPY Tax_kraken_to_phyloseq.R /mnt/
COPY Best_hits_search.R /mnt/
COPY Bin_summary.R /mnt/
COPY Bin_checkm_general_plot.R /mnt/
COPY Blobplot.R /mnt/
COPY Read_arg_norm.R /mnt/
COPY Contig_arg_unify.R /mnt/
COPY Tax_sourmash_to_phyloseq.R /mnt/
COPY SM_unify_report.R /mnt/
COPY ARG_blob_plot.R /mnt/
COPY Contig_arg_unify_by_groups.R /mnt/
RUN apt-get update && apt-get -y install wget && apt-get -y install curl && apt-get -y install parallel
RUN export PATH=/mnt/:"$PATH"
RUN apt-get update && apt-get install -y awscli
