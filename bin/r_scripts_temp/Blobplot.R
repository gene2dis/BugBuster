library(tidyverse)
library(rcartocolor)
library(cowplot) 

assign("File_name",
       list.files(path = '.',
                  pattern = '*Blob_table*'))

sample_names <- c()

for (name in File_name) {
  sample_names <- c(sample_names, unlist(strsplit(name,"_Blob_table"))[1])
}

for (samples in unique(sample_names)){
  
  assign("report_unify",
         paste(samples,"blob_table_data", sep = "_"))
  
  assign(report_unify,
         list.files(path = '.',
                    pattern = paste0(samples,'_Blob_table*')) %>%
           read_tsv(skip = 10) %>%
           mutate(Sample = samples))
}

all_blob_reports <- lapply(ls(pattern = "*blob_table_data"),
                           function(df_name){get(df_name)}) %>%
  reduce(rbind)

new_names <- all_blob_reports %>%
  select(contig_name = `# name`,length,GC,bam0,phylum = `phylum.t.9%s`) %>%
  group_by(phylum) %>%
  summarise(length_sum = sum(length)) %>%
  arrange(-length_sum) 

num_top <- ifelse(nrow(new_names) < 9, nrow(new_names), 9)

data_to_plot <- all_blob_reports %>%
  select(contig_name = `# name`,length,GC,bam0,phylum = `phylum.t.9%s`) %>%
  left_join(new_names, by = "phylum") %>%
  mutate(phylum = ifelse(length_sum <= pull(new_names[num_top,]), "other", phylum)) 

fill_colors <- carto_pal(length(unique(data_to_plot$phylum)), "Pastel")

color_df <- data.frame(
  phylum = unique(data_to_plot$phylum),
  colour = fill_colors
) %>%
  mutate(colour = ifelse(phylum == "other","gray40",
                         ifelse(phylum == "no-hit","white",colour)))

phylum_percent <- data_to_plot %>%
  select(contig_name, phylum) %>%
  group_by(phylum) %>%
  count() %>%
  ungroup() %>%
  summarise(percent = round(n/sum(n)*100,digits = 2),
            phylum = phylum)

library(plyr)

tryCatch({
  my_first_blob_plot <- ggplot(data_to_plot, aes(x = GC,
                         y = log(bam0,base = 10),
                         size = length)) +
  geom_point(colour = "gray30",
             pch = 21,
             alpha = 0.8,
             aes(fill = factor(phylum))) +
  scale_size(range = c(.1, 24), name="Contig size",
             breaks = c(round_any(max(data_to_plot$length)/12,10),
                        round_any(max(data_to_plot$length)/8,10),
                        round_any(max(data_to_plot$length)/4,10)),
             labels = c(paste0(round_any(max(data_to_plot$length)/12,10)," nt"),
                        paste0(round_any(max(data_to_plot$length)/8,10)," nt"),
                        paste0(round_any(max(data_to_plot$length)/4,10)," nt"))) +
  scale_fill_manual(values = setNames(color_df$colour, color_df$phylum),
                    guide = guide_legend(ncol=1)) +
  labs(fill = "Phylum",
       x = "GC content",
       y = "Coverage (log)")+
  guides(fill = guide_legend(override.aes = list(size=10))) +
  theme_bw(base_size=14) +
  scale_y_continuous(labels = c("",expression(paste("10"^"1")),
                                expression(paste("10"^"2")),
                                expression(paste("10"^"3")),
                                expression(paste("10"^"4"))),
                     breaks = c(0,1,2,3,4),
                     limits = c(0,4),
                     expand = c(0,0)) +
  scale_x_continuous(labels = c("0","20%","40%","60%","80%","100%"),
                     breaks = c(0,0.2,0.4,0.6,0.8,1),
                     limits = c(0,1),
                     expand = c(0,0)) +
  theme(panel.grid.major.x = element_line(color = "gray80"),
        panel.grid.major.y = element_line(color = "gray80"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position="right",
        legend.justification = "top",
        strip.background = element_rect(color="#CDDEFF", fill="#CDDEFF"),
        text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box = "vertical", 
        legend.direction = "vertical")
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  phylum_classified <- ggplot(phylum_percent) +
  geom_bar(aes(y = percent,
               x = phylum,
               fill = phylum),
           stat="identity",
           colour="black")  +
  geom_text(aes(y = percent,
                x = phylum,
                label = paste0(percent,"%")),
            vjust=-0.5,
            size = 6) +
  scale_fill_manual(values = setNames(color_df$colour, color_df$phylum),
                    guide = guide_legend(ncol=1)) +
  scale_y_continuous(labels = c("0","20%","40%","60%","80%","100%"),
                     breaks = c(0,20,40,60,80,100),
                     limits = c(0,100),
                     expand = c(0,0)) +
  theme_bw(base_size=14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray60",
                                          linetype = 2),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.line = element_line(color = "black"),
        legend.position= "none",
        legend.justification = NULL,
        strip.background = element_rect(color="#CDDEFF", fill="#CDDEFF"),
        text = element_text(size = 25),
        axis.text.x = element_text(size = 20,angle = 40,vjust = 1,hjust = 1),
        panel.ontop = TRUE) +
  labs(fill = "Phylum",
       x = NULL,
       y = NULL)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
phylum_to_cow <- plot_grid(phylum_classified, ggplot() + theme_minimal(),
          rel_widths = c(0.828, 0.172))
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
final_blob_plot <- plot_grid(my_first_blob_plot, phylum_to_cow,
          ncol = 1,
          rel_widths = c(1, 1),
          rel_heights = c(0.6,0.4))
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ggsave("Phylum_blob_plot.png", 
       plot = final_blob_plot, 
       dpi = 300,
       width = 20,
       height = 14)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

saveRDS(file = "Phylum_final_blobplot.RDS", object = final_blob_plot)
