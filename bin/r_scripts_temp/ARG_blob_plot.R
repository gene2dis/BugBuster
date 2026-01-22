library(tidyverse)
library(rcartocolor)
library(cowplot)

contig_arg_data <- read_csv("Contig_tax_and_arg_prediction.tsv") 

if(all(is.na(unique(contig_arg_data$ARG_gene)))){
  
  all_na <- TRUE
  
  contig_arg_data <- contig_arg_data %>%
    mutate(ARG_gene = ifelse(is.na(ARG_gene), "None",ARG_gene),
           predicted_ARG_class = ifelse(is.na(ARG_gene), "None",ARG_gene))
} else {
  
  all_na <- FALSE

}

new_names <- contig_arg_data %>%
  select(contig_name, contig_length, GC, coverage, predicted_ARG_class) %>%
  group_by(predicted_ARG_class) %>%
  summarise(length_sum = sum(contig_length)) %>%
  arrange(-length_sum)

if(length(new_names$predicted_ARG_class) < 9){
  max_num_pull <- length(new_names$predicted_ARG_class)
} else {
  max_num_pull <- 9
}

if(all_na){
  data_to_plot <- contig_arg_data %>%
    select(contig_name, contig_length, GC, coverage, predicted_ARG_class) %>%
    left_join(new_names, by = "predicted_ARG_class") %>%
    mutate(predicted_ARG_class = str_to_title(predicted_ARG_class)) %>%
    mutate(coverage = log(coverage,base = 10))
} else {
  data_to_plot <- contig_arg_data %>%
    select(contig_name, contig_length, GC, coverage, predicted_ARG_class) %>%
    left_join(new_names, by = "predicted_ARG_class") %>%
    mutate(predicted_ARG_class = ifelse(predicted_ARG_class %in%
                                          new_names$predicted_ARG_class[1:max_num_pull],
           predicted_ARG_class, "Other")) %>%
    mutate(predicted_ARG_class = str_to_title(predicted_ARG_class)) %>%
    filter(!is.na(predicted_ARG_class)) %>%
    mutate(coverage = log(coverage,base = 10)) %>%
    mutate(coverage = ifelse(coverage == -Inf, sort(unique(coverage))[2],coverage))
}



if(all_na){
  fill_colors <- "gray70"
} else {
  fill_colors <- carto_pal(length(unique(data_to_plot$predicted_ARG_class)), "Safe")
}


color_df <- data.frame(
  predicted_ARG_class = unique(data_to_plot$predicted_ARG_class),
  colour = fill_colors
)

predicted_ARG_class_percent <- data_to_plot %>%
  select(contig_name, predicted_ARG_class) %>%
  group_by(predicted_ARG_class) %>%
  count() %>%
  ungroup() %>%
  summarise(percent = round(n/sum(n)*100,digits = 2),
            predicted_ARG_class = predicted_ARG_class)

library(plyr)

tryCatch({

  my_first_blob_plot <- ggplot(data_to_plot, aes(x = GC,
                                                 y = coverage,
                                                 size = contig_length)) +
    geom_point(colour = "gray30",
               pch = 21,
               alpha = 0.8,
               aes(fill = factor(predicted_ARG_class))) +
    scale_size(range = c(.1, 24), name="Contig size",
               breaks = c(round_any(max(data_to_plot$contig_length)/12,10),
                          round_any(max(data_to_plot$contig_length)/8,10),
                          round_any(max(data_to_plot$contig_length)/4,10)),
               labels = c(paste0(round_any(max(data_to_plot$contig_length)/12,10)," nt"),
                          paste0(round_any(max(data_to_plot$contig_length)/8,10)," nt"),
                          paste0(round_any(max(data_to_plot$contig_length)/4,10)," nt"))) +
    scale_fill_manual(values = setNames(color_df$colour, color_df$predicted_ARG_class),
                      guide = guide_legend(ncol=1)) +
    labs(fill = "Predicted ARG Class",
         x = "GC content",
         y = "Coverage (log)")+
    guides(fill = guide_legend(override.aes = list(size=10))) +
    theme_bw(base_size=14) +
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
  ARG_classified <- ggplot(predicted_ARG_class_percent) +
    geom_bar(aes(y = percent,
                 x = predicted_ARG_class,
                 fill = predicted_ARG_class),
             stat="identity",
             colour="black")  +
    geom_text(aes(y = percent,
                  x = predicted_ARG_class,
                  label = paste0(percent,"%")),
              vjust=-0.5,
              size = 6) +
    scale_fill_manual(values = setNames(color_df$colour, color_df$predicted_ARG_class),
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
    labs(fill = "Predicted ARG Class",
         x = NULL,
         y = NULL)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ARG_to_cow <- plot_grid(ARG_classified, ggplot() + theme_minimal(),
                          rel_widths = c(0.828, 0.172))
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ARG_final_blob_plot <- plot_grid(my_first_blob_plot, ARG_to_cow,
                                   ncol = 1,
                                   rel_widths = c(1, 1),
                                   rel_heights = c(0.6,0.4))
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

tryCatch({
  ggsave("ARG_blob_plot.png",
         plot = ARG_final_blob_plot,
         dpi = 300,
         width = 20,
         height = 14)
}, error = function(e) {
  cat("Se ha producido un error:", conditionMessage(e), "\n")
})

saveRDS(file = "ARG_final_blobplot.RDS", object = ARG_final_blob_plot)
