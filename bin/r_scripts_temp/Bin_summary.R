library(tidyverse)

bin_depth <- list.files(path='.',pattern = "*bin_depth.tsv") %>%
  lapply(read_tsv) %>%
  bind_rows %>%
  mutate(sample = paste0(sample,"_avgcov")) %>%
  pivot_wider(names_from = sample,values_from = average_cov) 

bin_tmp <- bin_depth %>%
  select(bin_id, Total_contigs)

bin_quality <- list.files(path='.',pattern = "*quality_report.tsv") %>%
  read_tsv() %>%
  select(Name,Completeness,Contamination,Contig_N50,Genome_Size,Total_Coding_Sequences)

bin_tax <- list.files(path = '.',pattern = "*gtdbtk*") %>%
  read_tsv() %>%
  mutate(Domain = sapply(strsplit(classification, ";"),`[`,1),
         Phylum = sapply(strsplit(classification, ";"),`[`,2),
         Class = sapply(strsplit(classification, ";"),`[`,3),
         Order = sapply(strsplit(classification, ";"),`[`,4),
         Family = sapply(strsplit(classification, ";"),`[`,5),
         Genus = sapply(strsplit(classification, ";"),`[`,6),
         Species = sapply(strsplit(classification, ";"),`[`,7)) %>%
  mutate(Domain = sapply(strsplit(Domain, "__"),`[`,2),
         Phylum = sapply(strsplit(Phylum, "__"),`[`,2),
         Class = sapply(strsplit(Class, "__"),`[`,2),
         Order = sapply(strsplit(Order, "__"),`[`,2),
         Family = sapply(strsplit(Family, "__"),`[`,2),
         Genus = sapply(strsplit(Genus, "__"),`[`,2),
         Species = sapply(strsplit(Species, "__"),`[`,2)) %>%
  select(user_genome, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(Species = ifelse(is.na(Species), paste0(Genus," sp."),Species))

Bin_reports <- left_join(bin_quality, bin_tmp, by = join_by(Name == bin_id)) %>%
  left_join(bin_tax, by = join_by(Name == user_genome)) %>%
  left_join(bin_depth %>% select(!Total_contigs), by = join_by(Name == bin_id)) 

write.csv(Bin_reports, "Bin_summary.csv", quote = F,row.names = F)

