library(tidyverse)

assign("File_name",
       list.files(path = '.',pattern = "*KARGA*"))
sample_names <- c()

for (name in File_name) {
  sample_names <- c(sample_names, unlist(strsplit(name,"_all_reads_"))[1])
}

for (samples in unique(sample_names)){
  
  assign("KARGA_report_unify",
         paste(samples,"KARGA_report", sep = "_"))
  
  assign(KARGA_report_unify,
         list.files(path = '.',
                    pattern = paste0(samples,"_all_reads_KARGA_mappedGenes.csv")) %>%
           read_csv() %>%
           mutate(sample = samples))
  
  assign("KARGVA_report_unify",
         paste(samples,"KARGVA_report", sep = "_"))
  
  assign(KARGVA_report_unify,
         list.files(path = '.',
                    pattern = paste0(samples,"_all_reads_KARGVA_mappedGenes.csv")) %>%
           read_csv() %>%
           mutate(sample = samples))

  assign("ARGs_OAP_cells",
         paste(samples,"ARGs_OAP_cells_report",sep = "_"))

  sample_dirs <- list.dirs()
  sample_dirs <- sample_dirs[ grepl(paste0(samples,"_args_oap_s1_out"), sample_dirs) ]

  assign(ARGs_OAP_cells,
         list.files(path = sample_dirs,
                    pattern = "metadata.txt",
                    full.names = T) %>%
           read_tsv() %>%
           mutate(sample = sapply(strsplit(as.character(sample), split = "_tmp_"), `[`,1)) %>%
           group_by(sample) %>%
           summarise(nreads = sum(nRead),
                     n16S = sum(n16S),
                     nCell = sum(nCell)))  
}

All_KARGA_reports <- lapply(ls(pattern = "*_KARGA_report"),
                            function(df_name){get(df_name)}) %>%
  reduce(rbind)

All_KARGVA_reports <- lapply(ls(pattern = "*_KARGVA_report"),
                             function(df_name){get(df_name)}) %>%
  reduce(rbind)

All_ARGs_OAP_cells_reports <- lapply(ls(pattern = "*ARGs_OAP_cells_report"),
                             function(df_name){get(df_name)}) %>%
  reduce(rbind)

if (all(is.na(All_KARGA_reports$GeneIdx))) {
  KARGA_ARG_pred_data_filt <- All_KARGA_reports 
} else {
  KARGA_ARG_pred_data_filt <- All_KARGA_reports %>%
    mutate(Compound = sapply(strsplit(GeneIdx, split = "\\|"), `[`,2),
           Meg_drug = sapply(strsplit(GeneIdx, split = "\\|"), `[`,3),
           Other_Drug = sapply(strsplit(GeneIdx, split = "\\|"), `[`,4),
           AverageKMerDepth = as.numeric(AverageKMerDepth),
           Gen = sapply(strsplit(GeneIdx, split = "\\|"), `[`,5),
           PercentGeneCovered = as.numeric(sapply(strsplit(PercentGeneCovered, split = "%"), `[`,1))) %>%
    filter(PercentGeneCovered >= 90) %>%
    mutate(Meg_drug = ifelse(Compound == "Metals", "Metal resistance", Meg_drug)) %>%
    mutate(Meg_drug = ifelse(Compound == "Biocides", "Biocide resistance", Meg_drug)) %>%
    mutate(Meg_drug = ifelse(Compound == "Multi-compound", "Multi compound resistance", Meg_drug)) %>%
    mutate(Meg_drug = sapply(Meg_drug, function(palabra) {
      palabra <- trimws(palabra)
      paste(toupper(substring(palabra, 1, 1)),
            tolower(substring(palabra, 2)),
            sep = "")})) %>%
    select(sample,GeneIdx,Meg_drug,Compound,Gen,PercentGeneCovered,AverageKMerDepth) %>%
    mutate(Meg_drug = str_replace_all(Meg_drug, "_", " "))
}

if (all(is.na(All_KARGVA_reports$GeneIdx))) {
  KARGVA_ARG_pred_data_filt <- All_KARGVA_reports
} else {
  KARGVA_ARG_pred_data_filt <- All_KARGVA_reports %>%
    mutate(Mutation = sapply(strsplit(GeneIdx, split = "\\|"), `[`,2),
           Meg_drug = sapply(strsplit(GeneIdx, split = "\\|"), `[`,7),
           Drug = sapply(strsplit(GeneIdx, split = "\\|"), `[`, 6),
           Other_Drug = sapply(strsplit(GeneIdx, split = "\\|"), `[`,15),
           Other_Drug_2 = sapply(strsplit(GeneIdx, split = "\\|"), `[`,8),
           Organism_1 = sapply(strsplit(GeneIdx, split = "\\|"), `[`,4),
           Organism_2 = sapply(strsplit(GeneIdx, split = "\\|"), `[`,17),
           KmerSNPHits = as.numeric(sapply(strsplit(KmerSNPHits, split = "/"), `[`,1)),
           PercentGeneCovered = as.numeric(sapply(strsplit(PercentGeneCovered, split = "%"), `[`,1)),
           Gen_1 = sapply(strsplit(GeneIdx, split = "\\|"), `[`,9),
           Gen_2 = sapply(strsplit(GeneIdx, split = "\\|"), `[`,14),
           Gen_2 = sapply(strsplit(Gen_2, split = "_"),`[`,1),
           Gen_3 = sapply(strsplit(GeneIdx, split = "\\|"), `[`,4),
           Gen_3 = sapply(strsplit(Gen_3, split = " "),`[`,3),
           Gen_fusion = ifelse(!is.na(Gen_3), Gen_3,
                               ifelse(!is.na(Gen_2), Gen_2, Gen_1))) %>%
    filter(KmerSNPHits >= 2) %>%
    filter(PercentGeneCovered >= 80) %>%
    mutate(Organism_1 = sapply(strsplit(Organism_1, split = "resistance to"), `[`,2),
           Organism_2 = sapply(strsplit(Organism_2, split = "_resistant"), `[`,1)) %>%
    mutate(All_drug = ifelse(!is.na(Organism_1), Organism_1, 
                             ifelse(Other_Drug == "NA", Drug, Other_Drug))) %>%
    mutate(All_drug = ifelse(All_drug == "Multi-compound", "MULTIDRUG", All_drug),
           All_drug = ifelse(All_drug == "Drugs","MULTIDRUG", All_drug),
           All_drug = ifelse(All_drug == " rifampicin","RIFAMYCIN", All_drug),
           All_drug = ifelse(All_drug == " fusidic acid","FUSIDIC_ACID", All_drug),
           All_drug = ifelse(All_drug == " sulfonamides","SULFONAMIDE", All_drug)) %>%
    select(sample,GeneIdx,Drug_resistance = All_drug,Gen = Gen_fusion,PercentGeneCovered,AverageKMerDepth) %>%
    mutate(Drug_resistance = sapply(Drug_resistance, function(palabra) {
      palabra <- trimws(palabra)
      paste(toupper(substring(palabra, 1, 1)),
            tolower(substring(palabra, 2)),
            sep = "")})) 
}

reads_and_cells_report <- All_ARGs_OAP_cells_reports %>%
  select(sample,nRead = nreads,n16S,nCell)

KARGA_with_reads <- left_join(reads_and_cells_report, KARGA_ARG_pred_data_filt, by = "sample")
KARGVA_with_reads <- left_join(reads_and_cells_report, KARGVA_ARG_pred_data_filt, by = "sample")

KARGA_with_reads[] <- lapply(KARGA_with_reads, function(x) {
  if (is.list(x)) {
    x <- unlist(x)
    x <- as.character(NA)
  }
  return(x)
})

KARGVA_with_reads[] <- lapply(KARGVA_with_reads, function(x) {
  if (is.list(x)) {
    x <- unlist(x)
    x <- as.character(NA)
  }
  return(x)
})

if (all(is.na(All_KARGA_reports$GeneIdx))) {
  KARGA_norm <- KARGA_with_reads
} else {
  KARGA_norm <- KARGA_with_reads %>%
    group_by(sample,Gen) %>%
    mutate(CPM_ARGs = AverageKMerDepth*(1/(nRead/10^6))) %>%
    mutate(copies_per_cell = AverageKMerDepth/nCell)
}

if (all(is.na(All_KARGVA_reports$GeneIdx))) {
  KARGVA_norm <- KARGVA_with_reads 
} else {
  KARGVA_norm <- KARGVA_with_reads %>%
    group_by(sample,Gen) %>%
    mutate(CPM_ARGs = AverageKMerDepth*(1/(nRead/10^6))) %>%
    mutate(copies_per_cell = AverageKMerDepth/nCell)
}

write.csv(KARGA_norm, file = "KARGA_norm.csv", row.names = F, quote = F)
write.csv(KARGVA_norm, file = "KARGVA_norm.csv", row.names = F, quote = F)
