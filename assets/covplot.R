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

# written by Carl Suster @arcesu & Winkie Fong @lilwinx

# input data
args <- commandArgs(trailingOnly = TRUE)
cov_file <- args[1]
cov <- read_tsv(cov_file)
outdir <- args[2]
#outdir <- "/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/coverage_ref.test"
#filter_coords <- read_tsv("/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/coverage_contigs_filter_coords.txt", skip = 3)
#cov <- read_tsv("/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/coverage_ref.depth.txt")
outname <- "coverageplot.png"

# filter the filter_coords and cov file if it contains more than one unique value
split_data <- split(cov, cov$ID)
cov_var_names <- paste0("cov_", names(split_data))
for (i in seq_along(split_data)) {
  assign(cov_var_names[i], split_data[[i]])
}

grouped_data_list <- list()
for (entry in unique(unlist(lapply(split_data, function(df) df$ID)))) {
  data_names <- c(paste0("cov_", entry))
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
plot_binned_coverage <- function(coverage_df, ref_name) {
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
    scale_x_continuous(labels = comma) +
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


# creating the new list to loop through based on entries.
for (entry_data_names in grouped_data_list) {
  entry_name <- gsub("cov_", "", entry_data_names[1])
  cov_name <- get(entry_data_names[1])
  cov_df <- create_binned_coverage(cov_name)
  
  plot <- plot_binned_coverage(cov_df, entry_name)
  
  cov_df_name <- paste0("cov_df_", entry_name)
  cov_plot_name <- paste0("cov_plot_", entry_name)
  combined_plot_name <- paste0("combined_plot_", entry_name)
  
  assign(cov_df_name, cov_df)
  assign(cov_plot_name, plot)
  assign(combined_plot_name, plot)
  
  ggsave(paste(outdir, entry_name, outname, sep = "_"),
         width = 1500,
         height = 500,
         units = "px",
         dpi = 300,
         scale = 2,
         bg = "transparent")
  print("Output file located in:")
  print(paste(outdir, entry_name, outname, sep = ""))
}
