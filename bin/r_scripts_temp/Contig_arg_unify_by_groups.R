library(tidyverse)

assign("File_name",
       list.files(path = '.',
                  pattern = '*Blob_table*'))

sample_names <- c()

for (name in File_name) {
  sample_names <- c(sample_names, unlist(strsplit(name,"_Blob_table"))[1])
}

for (samples in unique(sample_names)){
  
  sample_group <- unlist(strsplit(samples,"_from_"))[1]
  sample_name <- unlist(strsplit(samples,"_from_"))[2]
  
  assign("report_unify",
         paste(samples,"blob_table_data", sep = "_"))
  
  assign(report_unify,
         list.files(path = '.',
                    pattern = paste0(samples,'_Blob_table*')) %>%
           read_tsv(skip = 10) %>%
           mutate(Sample_group = sample_group,
                  Sample_id = sample_name))
}

assign("File_name_arg",
       list.files(path = '.',
                  pattern = '*_contigs_deep_arg.out.mapping.ARG'))

sample_names_for_deeparg <- c()

for (name in File_name_arg) {
  sample_names_for_deeparg <- c(sample_names_for_deeparg, unlist(strsplit(name,"_contigs_deep_arg.out.mapping.ARG"))[1])
}

for (samples in unique(sample_names_for_deeparg)){
  
  assign("deep_arg_data",
         paste(samples,"deep_arg_data", sep = "_"))
  
  assign(deep_arg_data,
         list.files(path = '.',
                    pattern = paste0(samples,'_contigs_deep_arg.out.mapping.ARG')) %>%
         read_tsv() %>%
         mutate(Sample_group = samples))
  
}

all_blob_reports <- lapply(ls(pattern = "*blob_table_data"),
                           function(df_name){get(df_name)}) %>%
  reduce(rbind) %>%
  select(Sample_group, Sample_id, contig_name = `# name`,
         contig_length = length, GC,
         coverage = bam0,
         superkingdom = `superkingdom.t.6%s`,
         phylum = `phylum.t.9%s`,
         order = `order.t.12%s`,
         family = `family.t.15%s`,
         genus = `genus.t.18%s`,
         species = `species.t.21%s`) %>%
  mutate(Sample_id = paste0(Sample_id,"_avg_cov")) 

all_arg_reports <- lapply(ls(pattern = "*_deep_arg_data"),
                          function(df_name){get(df_name)}) %>%
  reduce(rbind) %>%
  mutate(new_read_id_1 = sapply(strsplit(read_id, "_"), `[`,1),
         new_read_id_2 = sapply(strsplit(read_id, "_"), `[`,2),
         contig_name = paste(new_read_id_1,new_read_id_2,sep = "_")) %>%
  select(Sample_group, contig_name, ARG_gene = `#ARG`,predicted_ARG_class = `predicted_ARG-class`,
         best_hit = `best-hit`, probability, identity,alignment_evalue = `alignment-evalue` ,counts)

all_data <- left_join(all_blob_reports,
                      all_arg_reports,
                      by = c("Sample_group","contig_name")) %>%
  select(Sample_group,Sample_id,contig_name,contig_length,GC,coverage,
         superkingdom,phylum,order,family,genus,species,
         ARG_gene,predicted_ARG_class,
         best_hit,probability,identity,alignment_evalue,counts)

write.csv(all_data,
          "Contig_tax_and_arg_prediction.csv",
          quote = F,
          row.names = F)
