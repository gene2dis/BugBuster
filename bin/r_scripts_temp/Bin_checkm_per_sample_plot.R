library(tidyverse)

assign("File_name",
       list.files(path = '.'))
sample_names <- c()

for (name in File_name) {
  sample_names <- c(sample_names, unlist(strsplit(name,"_*_quality_report"))[1])
  print(unlist(strsplit(name,"_"))[1])
}

for (samples in unique(sample_names)){
  for (binner in c("metabat","semibin","comebin","metawrap")) {
    
    assign("report_unify",
           paste(samples,binner,"quality_report", sep = "_"))
    
    assign(report_unify,
           list.files(path = '.',
                      pattern = paste0(samples,".+?",binner,".+")) %>%
             read_tsv() %>%
             mutate(Binner = binner,
                    Sample = samples))
  }
}

all_quality_reports <- lapply(ls(pattern = "*_quality_report"),
                              function(df_name){get(df_name)}) %>%
  reduce(rbind) 
  
all_quality_reports_qf <- all_quality_reports %>%
  mutate(`Quality classify` = ifelse(Completeness >= 50 & Contamination <= 5,
                                     "zMid", "Low"),
         `Quality classify` = ifelse(Completeness >= 90 & Contamination <= 5,
                                     "zzHigh", `Quality classify`)) %>%
  select(Binner, `Quality classify`,Sample) %>%
  group_by(Sample,Binner) %>%
  count(`Quality classify`, name = "Bins count") %>%
  mutate(Binner = str_to_title(Binner)) %>%
  mutate(type_of_MAG = ifelse(Binner != "Metawrap", "Raw MAGs", "Refined MAGs"))

all_quality_reports_qf$Binner <- factor(all_quality_reports_qf$Binner,
                                        levels = c("Metabat","Comebin","Semibin","Metawrap"))

tryCatch({
  Quality_plot <- ggplot(all_quality_reports_qf, aes(x = Binner, y = `Bins count`, fill = `Quality classify`)) +
    geom_col(position = "dodge",
             colour = "black",
             linewidth = 1) +
    facet_wrap(~type_of_MAG) +
    scale_fill_manual(values=c("#d16678","#6baed6","#78b33e"),
                      labels = c("Low \n(Completeness < 50%, Contamination > 5%)",
                                 "Mid \n(Completeness => 50% < 90%, Contamination <= 5%)",
                                 "High \n(Completeness => 90%, Contamination <= 5%)")) +
    xlab("Binner") +
    scale_y_continuous(n.breaks = 12) +
    theme_bw(base_size = 14) +
    theme(axis.line = element_line(color = "black"),
          legend.justification = "top",
          legend.text = element_text(size=12),
          legend.spacing.y = unit(40, 'cm'),
          strip.background = element_rect(color="#CDDEFF", fill="#CDDEFF"),
          panel.border = element_blank(),
          text = element_text(size = 20),
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(1, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(color = "black",
                                     angle = 30, hjust = 1, vjust = 1, size = 18)) 
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ggsave("Sample_bins_quality_plot.png",
         plot = Quality_plot,
         dpi = 300,
         width=16,
         height=8)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

