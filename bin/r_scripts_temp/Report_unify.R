library("tidyverse")
library("plyr")

args = commandArgs(trailingOnly=TRUE)

if (args[1] != "none") {
	for (arg in args) {	
	  report_unify <- paste("bowtie", arg, sep = "_")
	  assign(report_unify,
	         list.files(path = '.',
	                    pattern = paste0("*",arg,"*")) %>%
	           map(read_tsv, col_types = "cd") %>%
	           bind_rows)
	}

	fastp_reports <- list.files(path='.',pattern = "*fastp_report.tsv") %>%
	  map(read_tsv, col_types = "cdd") %>%
	  bind_rows

	get_my_df <- function(df_name) {
	  get(df_name)
	}

	all_bowtie_reports <- lapply(ls(pattern = "*bowtie*"),get_my_df) %>%
	  reduce(dplyr::left_join, by = 'Id')

	all_reports <- left_join(fastp_reports, all_bowtie_reports) %>% mutate_all(~replace(., is.na(.), 0))

	without_id <- all_reports[,2:ncol(all_reports)]
	without_id_order <- without_id[,names(sort(colSums(without_id), decreasing = TRUE))]

	all_reports_order <- all_reports[,c("Id",names(sort(colSums(without_id), decreasing = TRUE)))]

        data_plot <- all_reports_order %>%
          pivot_longer(names_to = "Process",
                       values_to = "Reads",
                       cols = -c("Id")) %>%
          mutate(Read_type = ifelse(grepl(pattern = "singletons", Process),"Singleton","Paired")) %>%
          mutate(Process = str_replace(Process, " singletons", "")) %>%
          mutate(Process = ifelse(grepl(pattern = "Raw", Process), "Raw reads",Process)) %>%
          mutate(Process = ifelse(grepl(pattern = args[1], Process), paste0(args[1]," reads removal"),
                                  ifelse(grepl(pattern = args[2], Process), paste0(args[2]," reads removal"), Process)))

	data_plot$Process <- factor(data_plot$Process, levels = c("Raw reads","Fastp",paste0(args[1]," reads removal"),paste0(args[2]," reads removal")))

	max_round <- ceiling(max(data_plot$Reads) / 5000000)  *5000000

	tryCatch({
	boxplot_reads <- ggplot(data_plot, aes(x=Process, y=Reads, color=Process)) +
	    xlab("Process") +
	    ylab("Total reads") +
	    scale_y_continuous(labels = scales::number_format(scale = 1e-6, suffix = "M"),
	                       breaks = seq(0, max_round, 5000000),
	                       limits = c(0, max_round)) +
	    facet_wrap(~Read_type) +
	    geom_boxplot(outlier.colour="blue") +
	    geom_jitter(size=3, alpha=0.4) +
	    theme_bw() +
	    theme(axis.line = element_line(color = "black"),
	          legend.position = "none",
	          strip.background = element_rect(color="#CDDEFF", fill="#CDDEFF"),
	          text = element_text(size = 20),
	          axis.text.x = element_text(angle=60,hjust = 1,vjust = 1),
	          axis.ticks.length = unit(0.1,units = "cm"))
	}, error = function(e) {
	  cat("Se ha producido un error:", conditionMessage(e), "\n")
	})

	tryCatch({
	  ggsave("Box_plot_reads.png",
	         plot = boxplot_reads,
	         dpi = 300,
	         width=16,
	         height=8)
	}, error = function(e) {
	  cat("Se ha producido un error:", conditionMessage(e), "\n")
	})

	write.csv(all_reports_order,
	          "Reads_report.csv",
	          row.names = F,
	          quote = F)
} else {
        fastp_reports <- list.files(path='.',pattern = "*fastp_report.tsv") %>%
          map(read_tsv) %>%
          bind_rows

        write.csv(fastp_reports,
                  "Reads_report.csv",
                  row.names = F,
                  quote = F)
}
