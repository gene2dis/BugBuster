library(tidyverse)
library(phyloseq)
library(microViz)

db_used = commandArgs(trailingOnly=TRUE)

# read in sourmash gather & taxonomy results ------------------------------

# read in the sourmash taxonomy results from all samples into a single data frame
sourmash_taxonomy_results <- Sys.glob("*smgather*") %>%
  map_dfr(read_csv, col_types = "ddddddddcccddddcccdc") %>%
  mutate(name = gsub(" .*", "", name)) 

# We need two tables: a tax table and an "otu" table. 
# The tax table will hold the taxonomic lineages of each of our gather matches.
# To make this, we'll make a table with two columns: the genome match and the lineage of the genome.
# The "otu" table will have the counts of each genome in each sample.
# We'll call this our gather_table.

tax_table <- sourmash_taxonomy_results %>%
  select(name, lineage) %>%
  distinct() %>%
  separate(lineage, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  column_to_rownames("name") %>%
  mutate_at(.vars = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),.funs = function(x){gsub(".__", "", x)})

gather_table <- sourmash_taxonomy_results %>% 
  mutate(n_unique_kmers = (unique_intersect_bp / scaled) * average_abund) %>% # calculate the number of uniquely matched k-mers and multiply it by average k-mer abundance
  select(query_name, name, n_unique_kmers) %>% # select only the columns that have information we need
  pivot_wider(id_cols = name, names_from = query_name, values_from = n_unique_kmers) %>% # transform to wide format
  replace(is.na(.), 0) %>% # replace all NAs with 0 
  column_to_rownames("name") # move the metagenome sample name to a rowname

# create a phyloseq object ------------------------------------------------

physeq <- phyloseq(otu_table(gather_table, taxa_are_rows = T),
                   tax_table(as.matrix(tax_table)))

saveRDS(physeq,paste0("sourmash_",db_used,"_ps.RDS"))

tryCatch({
  assign("sourmash_phylum_bar", get("physeq") %>%
           comp_barplot(tax_level = "Phylum",
                        n_taxa = 10,
                        bar_width = 0.9,
                        bar_outline_colour = "black") +
           ggtitle(paste0("Relative abundance Phylum plot from ",toupper(db_used))) +
           scale_y_continuous(expand = c(0,0),
                              labels = scales::percent_format(scale=100)) +
           coord_flip() +
           ylab("Relative abundance") +
           theme(axis.line = element_line(color = "black"),
                 panel.grid = element_blank(),
                 legend.justification = "top",
                 legend.text = element_text(size=12),
                 legend.spacing.y = unit(40, 'cm'),
                 strip.background = element_rect(color="#CDDEFF", fill="#CDDEFF"),
                 text = element_text(size = 20),
                 legend.key.width = unit(1, "cm"),
                 legend.key.height = unit(1, "cm"),
                 legend.margin = margin(0, 0, 0, 0),
                 axis.text.y = element_text(size = 18),
                 axis.text.x = element_text(color = "black", size = 18)))
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ggsave(paste0("sourmash_",db_used,"_phylum_bar.png"),
         plot = get("sourmash_phylum_bar"),
         dpi = 300,
         width=16,
         height=8)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  assign("sourmash_fam_bar", get("physeq") %>%
           comp_barplot(tax_level = "Family",
                        n_taxa = 10,
                        bar_width = 0.9,
                        bar_outline_colour = "black") +
           ggtitle(paste0("Relative abundance Family plot from ",toupper(db_used))) +
           scale_y_continuous(expand = c(0,0),
                              labels = scales::percent_format(scale=100)) +
           coord_flip() +
           ylab("Relative abundance") +
           theme(axis.line = element_line(color = "black"),
                 panel.grid = element_blank(),
                 legend.justification = "top",
                 legend.text = element_text(size=12),
                 legend.spacing.y = unit(40, 'cm'),
                 strip.background = element_rect(color="#CDDEFF", fill="#CDDEFF"),
                 text = element_text(size = 20),
                 legend.key.width = unit(1, "cm"),
                 legend.key.height = unit(1, "cm"),
                 legend.margin = margin(0, 0, 0, 0),
                 axis.text.y = element_text(size = 18),
                 axis.text.x = element_text(color = "black", size = 18)))
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ggsave(paste0("sourmash_",db_used,"_fam_bar.png"),
         plot = get("sourmash_fam_bar"),
         dpi = 300,
         width=16,
         height=8)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  assign("sourmash_genus_bar", get("physeq") %>%
           comp_barplot(tax_level = "Genus",
                        n_taxa = 10,
                        bar_width = 0.9,
                        bar_outline_colour = "black") +
           ggtitle(paste0("Relative abundance Genus plot from ",toupper(db_used))) +
           scale_y_continuous(expand = c(0,0),
                              labels = scales::percent_format(scale=100)) +
           coord_flip() +
           ylab("Relative abundance") +
           theme(axis.line = element_line(color = "black"),
                 panel.grid = element_blank(),
                 legend.justification = "top",
                 legend.text = element_text(size=12),
                 legend.spacing.y = unit(40, 'cm'),
                 strip.background = element_rect(color="#CDDEFF", fill="#CDDEFF"),
                 text = element_text(size = 20),
                 legend.key.width = unit(1, "cm"),
                 legend.key.height = unit(1, "cm"),
                 legend.margin = margin(0, 0, 0, 0),
                 axis.text.y = element_text(size = 18),
                 axis.text.x = element_text(color = "black", size = 18)))
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ggsave(paste0("sourmash_",db_used,"_genus_bar.png"),
         plot = get("sourmash_genus_bar"),
         dpi = 300,
         width=16,
         height=8,
         device=png)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  assign("sourmash_species_bar", get("physeq") %>%
           comp_barplot(tax_level = "Species",
                        n_taxa = 10,
                        bar_width = 0.9,
                        bar_outline_colour = "black") +
           ggtitle(paste0("Relative abundance Genus plot from ",toupper(db_used))) +
           scale_y_continuous(expand = c(0,0),
                              labels = scales::percent_format(scale=100)) +
           coord_flip() +
           ylab("Relative abundance") +
           theme(axis.line = element_line(color = "black"),
                 panel.grid = element_blank(),
                 legend.justification = "top",
                 legend.text = element_text(size=12),
                 legend.spacing.y = unit(40, 'cm'),
                 strip.background = element_rect(color="#CDDEFF", fill="#CDDEFF"),
                 text = element_text(size = 20),
                 legend.key.width = unit(1, "cm"),
                 legend.key.height = unit(1, "cm"),
                 legend.margin = margin(0, 0, 0, 0),
                 axis.text.y = element_text(size = 18),
                 axis.text.x = element_text(color = "black", size = 18)))
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ggsave(paste0("sourmash_",db_used,"_species_bar.png"),
         plot = get("sourmash_species_bar"),
         dpi = 300,
         width=16,
         height=8,
         device=png)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})
