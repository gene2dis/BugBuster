library("tidyverse")
library("plyr")

args = commandArgs(trailingOnly=TRUE)

for (arg in args) {
  reports <- paste(arg,"tax_report",sep = "_")
  assign(reports, list.files(path='.',
                                     pattern = paste0("*",arg,"_report.tsv")) %>%
    map(read_tsv, col_types = "ccdd") %>%
    bind_rows)
}


get_my_df <- function(df_name) {
  get(df_name)
}

reads_report <- read_csv("Reads_report.csv")

if (length(args) > 1) {
  all_tax_reports <- lapply(ls(pattern = "*tax_report*"),get_my_df) %>%
    bind_rows

  all_reports <- left_join(reads_report, all_tax_reports)
} else {
  all_reports <- left_join(reads_report, get(reports))
}

if (length(args) > 1) {
  tax_data <- all_tax_reports %>%
    pivot_longer(-c(Id,colnames(all_tax_reports)[2]), names_to = "Tax results",values_to = "Percent")
} else {
  tax_data <- get(reports) %>%
    pivot_longer(-c(Id,colnames(get(reports))[2]), names_to = "Tax results",values_to = "Percent")
}

tax_data_plot <- tax_data

tryCatch({
  tax_plot <- ggplot(tax_data_plot, aes(fill = `Tax results`, x = Percent, y = Id)) +
    geom_bar(position="stack", stat="identity", colour="black") +
    ylab("Samples") +
    theme_bw() +
    scale_fill_manual(values=c('#6baed6', '#d16678')) +
    scale_x_continuous(expand = c(0,0),
                       labels = scales::percent_format(scale=100)) +
    xlab("Reads Classified") +
    theme(axis.line = element_line(color = "black"),
          panel.border = element_rect(color = "black"),
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
          axis.text.x = element_text(color = "black", size = 18))

}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ggsave("sourmash_tax_classified_reads.png",
         plot = tax_plot,
         dpi = 300,
         width=16,
         height=8)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

write.csv(all_reports,
          "Reads_report.csv",
          row.names = F,
          quote = F)


