library(tidyverse)

assign("File_name",
       list.files(path = '.'))
sample_names <- c()

for (name in File_name) {
  sample_names <- c(sample_names, unlist(strsplit(name,"_*_quality_report"))[1])
}

for (samples in unique(sample_names)){
  
  assign("binner",
         unlist(strsplit(samples,".+_"))[2]) 
  
  assign("report_unify",
           paste(samples,"quality_report", sep = "_"))
  
  assign(report_unify,
           list.files(path = '.',
                      pattern = paste0(samples,".+")) %>%
             read_tsv() %>%
             mutate(Binner = binner,
                    Sample = sapply(strsplit(samples, paste0("_",binner)),`[`,1)))
}

all_quality_reports <- lapply(ls(pattern = "*_quality_report"),
                              function(df_name){get(df_name)}) %>%
  reduce(rbind)

all_quality_reports_qf <- all_quality_reports %>%
  mutate(`Quality classify` = ifelse(Completeness >= 50 & Contamination <= 5,
                                     "Mid", "Low"),
         `Quality classify` = ifelse(Completeness >= 90 & Contamination <= 5,
                                     "High", `Quality classify`)) %>%
  select(Binner, `Quality classify`,Completeness,Contamination,Sample,Name) %>%
  mutate(Binner = str_to_title(Binner)) %>%
  mutate(Bin_type = ifelse(Binner == "Metawrap", "Refined MAGs","Raw MAGs"))

all_quality_reports_qf$`Quality classify` <- factor(all_quality_reports_qf$`Quality classify`,
                                                    levels = c("Low","Mid","High"))

Mag_summary_table <- all_quality_reports_qf %>%
  group_by(Binner, `Quality classify`) %>%
  count(`Quality classify`, name = "Total MAGs")

All_data <- all_quality_reports_qf %>%
  select(-Bin_type) %>%
  select(Sample, Binner,`Bin name` = Name,`Quality classify`,Completeness,Contamination) 

refined_data <- all_quality_reports_qf %>%
  select(-Bin_type) %>%
  select(Sample, Binner,`Bin name` = Name,`Quality classify`,Completeness,Contamination) %>%
  filter(Binner == "Metawrap") 

tryCatch({
  General_mag_plot <- ggplot(all_quality_reports_qf, aes(x = Contamination,
                                                         y = Completeness,
                                                         colour = `Quality classify`)) +
    geom_point(alpha=0.5,
               size=4) +
    geom_vline(xintercept = 5,
               linetype="dashed",
               colour="black",
               size=0.8) +
    geom_hline(yintercept = 90,
               linetype="dashed",
               colour="black",
               size=0.8) +
    geom_hline(yintercept = 50,
               linetype="dashed",
               colour="black",
               size=0.8) +
    facet_wrap(~ Bin_type) +
    labs(colour = "Quality classification") +
    scale_colour_manual(values=c("#d16678","#6baed6","#78b33e"),
                        labels = c("Low \n(Completeness < 50%, Contamination > 5%)",
                                   "Mid \n(Completeness => 50% < 90%, Contamination <= 5%)",
                                   "High \n(Completeness => 90%, Contamination <= 5%)")) +
    scale_y_continuous(n.breaks = 12,
                       limits = c(0,100),
                       labels = scales::percent_format(scale = 1)) +
    scale_x_continuous(limits = c(0,20),
                       labels =scales::percent_format(scale = 1)) +
    theme_bw(base_size = 14) +
    theme(axis.line = element_line(color = "black"),
          legend.justification = "top",
          legend.text = element_text(size=18),
          legend.spacing.y = unit(40, 'cm'),
          strip.background = element_rect(color="#CDDEFF", fill="#CDDEFF"),
          text = element_text(size = 20),
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(1.5, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(color = "black",size = 18),
          legend.key.size = unit(1.5, 'cm'))
  
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ggsave("Total_bins_quality_plot.png",
         plot = General_mag_plot,
         dpi = 300,
         width=16,
         height=8)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

write.csv(Mag_summary_table, "Mag_quality_summary.csv",
          quote = F, row.names = F)

write.csv(All_data, "All_Mag_quality_table.csv",
          quote = F, row.names = F)

write.csv(refined_data, "Refined_Mag_quality_table.csv",
          quote = F, row.names = F)
