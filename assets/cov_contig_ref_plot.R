list.of.packages <- c(
  "ggplot2",
  "tidyverse",
  "ivs", 
  "patchwork", 
  "purrr", 
  "grid",
  "dplyr",
  "scales"
  )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(tidyverse)
library(ivs)
library(patchwork)
library(purrr)
library(grid)
library(scales)
library(dplyr)

# written by Carl Suster @arcesu

# input data
args <- commandArgs(trailingOnly = TRUE)
filter_coords_file <- args[1]
filter_coords <- read_tsv(filter_coords_file, skip = 3)
cov_file <- args[2]
cov <- read_tsv(cov_file)
outdir <- args[3]
#filter_coords <- read_tsv("/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/coverage_contigs_filter_coords.txt", skip = 3)
#cov <- read_tsv("/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/coverage_ref.depth.txt")
outname <- "_coverageplot.png"

# clean up contigs file to collapse overlapping regions
merge_overlap <- function(fc_df) {
  endpoint <- fc_df %>%
    rename_with(.fn = ~ stringr::str_remove_all(.x, "^\\[|\\]$")) %>%
    filter(`% IDY` > 87) %>%
    select(S1, E1) %>%
    arrange(S1)
  
  merged <- iv(endpoint$S1, endpoint$E1 + 1) %>% iv_groups()
  merged_df <- tibble(start = iv_start(merged), end = iv_end(merged) - 1)
  
  return(merged_df)
}

# filter the filter_coords and cov file if it contains more than one unique value
split_data <- split(cov, cov$ID)
cov_var_names <- paste0("cov_", names(split_data))
for (i in seq_along(split_data)) {
  assign(cov_var_names[i], split_data[[i]])
}

filter_coords$`[TAGS]` <- sapply(strsplit(filter_coords$`[TAGS]`, "\t"), function(x) x[1])
coord_split_data <- split(filter_coords, filter_coords$`[TAGS]`)
coords_var_names <- paste0("fc_", names(coord_split_data))
for (i in seq_along(coord_split_data)) {
  assign(coords_var_names[i], coord_split_data[[i]])
}

grouped_data_list <- list()
for (entry in unique(unlist(lapply(split_data, function(df) df$ID)))) {
  data_names <- c(paste0("fc_", entry), paste0("cov_", entry))
  grouped_data_list[[entry]] <- data_names
}

# commentable section to bin coverage figure
create_binned_coverage <- function(bin_df) {
  binwidth <- nrow(bin_df) / 10000
  bin_starts <- seq(1L, nrow(bin_df) - binwidth, binwidth)
  
  binned_cov <- tibble(
    Position = bin_starts + binwidth / 2,
    Coverage = map2_dbl(bin_starts, bin_starts + binwidth, ~ median(bin_df$Coverage[.x:.y]))
  )
  
  return(binned_cov)
}

# draw coverage plot
plot_binned_coverage <- function(coverage_df, filter_coords_df, ref_name) {
  p <- ggplot(
    coverage_df, 
    #cov, # switch to this for no bin!
    aes(x = Position, y = Coverage)) +
    geom_col(fill = "#5D3A9B") +
    labs(
      title = ref_name,
      x = NULL,
      y = "Fold10 Coverage"
    ) +
    theme_classic() +
    scale_x_continuous(labels = comma, limits = c(1, max(coverage_df$Position, filter_coords_df$end))) +
    scale_y_continuous(trans = 'log10') + 
    theme(
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 8),
      panel.background = element_rect(fill = "transparent", colour = NA), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA)
    )
  
  return(p)
}

# draw contigs plot
plot_contigs <- function(fc_df, cov_df, ref_name) {
  p <- ggplot(fc_df) +
    geom_segment(aes(x = start, xend = end, y = 0, yend = 0), colour = "#E66100", linewidth = 10) +
    annotation_custom(grid::textGrob("Contigs", x = grid::unit(0, "mm"), hjust = 0, gp = grid::gpar(fontsize = 8))) +
    labs(
      x = "Base position in genome",
      y = NULL
    ) +
    scale_x_continuous(labels = NULL, limits = c(1, max(cov_df$Position, fc_df$end))) +
    scale_y_continuous(labels = NULL) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA_character_),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA_character_),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA)
    )
  
  return(p)
}

# creating the new list to loop through based on entries.
for (entry_data_names in grouped_data_list) {
  entry_name <- gsub("fc_", "", entry_data_names[1])
  fc_name <- get(entry_data_names[1])
  cov_name <- get(entry_data_names[2])
  fc_df <- merge_overlap(fc_name)
  cov_df <- create_binned_coverage(cov_name)
  
  cov_plot <- plot_binned_coverage(cov_df, fc_df, entry_name)
  fc_plot <- plot_contigs(fc_df, cov_df, entry_name)
  plot <- cov_plot / fc_plot + 
      plot_layout(ncol = 1, heights = c(10, 1)) + 
      plot_annotation(theme = theme(
        plot.background = element_rect(fill = "transparent", color = NA)
        ))
  
  
  fc_df_name <- paste0("fc_df_", entry_name)
  cov_df_name <- paste0("cov_df_", entry_name)
  cov_plot_name <- paste0("cov_plot_", entry_name)
  fc_plot_name <- paste0("fc_plot_", entry_name)
  combined_plot_name <- paste0("combined_plot_", entry_name)
  
  assign(fc_df_name, fc_df)
  assign(cov_df_name, cov_df)
  assign(cov_plot_name, cov_plot)
  assign(fc_plot_name, fc_plot)
  assign(combined_plot_name, plot)
  
  ggsave(paste(outdir, entry_name, outname, sep = ""),
      width = 1500,
      height = 500,
      units = "px",
      dpi = 300,
      scale = 2,
      bg = "transparent")
  print("Output file located in:")
  print(paste(outdir, entry_name, outname, sep = ""))
}
