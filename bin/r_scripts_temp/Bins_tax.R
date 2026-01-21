library(tidyverse)
library(rcartocolor)

assign("File_name",
       list.files(path = '.',pattern = "*gtdbtk*"))
sample_names <- c()
for (name in File_name) {
  sample_names <- c(sample_names, unlist(strsplit(name,"_gtdbtk_"))[1])
}

for (samples in unique(sample_names)){

  assign("report_unify",
         paste(samples,"bin_tax_report", sep = "_"))

  assign(report_unify,
         list.files(path = '.',
                    pattern = paste0(samples,"_*")) %>%
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
    mutate(sample = samples) %>%
    select(sample, bin_name = `user_genome`, closest_genome_reference, Domain, Phylum, Class, Order, Family, Genus, Species))
}

all_tax_bin_reports <- lapply(ls(pattern = "*_bin_tax_report"),
                              function(df_name){get(df_name)}) %>%
  reduce(rbind) %>%
  select(bin_name,Phylum) %>%
  count(Phylum,name = "Bin count")

all_tax_bin_reports_csv <- lapply(ls(pattern = "*_bin_tax_report"),
                              function(df_name){get(df_name)}) %>%
  reduce(rbind)

write.csv(all_tax_bin_reports_csv, "MAGs_tax_summary.csv",
          quote = F, row.names = F)

fill_colors <- colorRampPalette(carto_pal(12,"Pastel"))(length(unique(all_tax_bin_reports$Phylum)))

color_df <- data.frame(
  phylum = unique(all_tax_bin_reports$Phylum),
  colour = fill_colors
)

tryCatch({
  bin_tax_plot <- ggplot(all_tax_bin_reports) +
    geom_bar(aes(y = `Bin count`,
                 x = Phylum,
                 fill = Phylum),
             stat="identity",
             colour="black")  +
    geom_text(aes(y = `Bin count`,
                  x = Phylum,
                  label = `Bin count`),
              vjust=-0.5,
              size = 6) +
    scale_fill_manual(values = setNames(color_df$colour, color_df$phylum),
                      guide = guide_legend(ncol=1)) +
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
         y = "MAGs Count")
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ggsave("MAGs_tax_plot.png",
         plot = bin_tax_plot,
         dpi = 300,
         width=16,
         height=8)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})
