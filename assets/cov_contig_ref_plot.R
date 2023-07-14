list.of.packages <- c(
  "ggplot2",
  "tidyverse",
  "ivs", 
  "patchwork", 
  "purrr", 
  "grid"
  )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(tidyverse)
library(ivs)
library(patchwork)
library(purrr)
library(grid)

# written by Carl Suster @arcesu

# input data
args <- commandArgs(trailingOnly = TRUE)
filter_coords_file <- args[1]
filter_coords <- read_tsv(filter_coords_file, skip = 3)
cov_file <- args[2]
cov <- read_tsv(cov_file)
outname <- args[3]
#filter_coords <- read_tsv("/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/coverage_contigs_filter_coords.txt", skip = 3)
#cov <- read_tsv("/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/coverage_ref.depth.txt")
#outname <- "output_coverageplot.png"

# clean up contigs file to collapse overlapping regions
endpoint <- filter_coords %>%
  rename_with(.fn = ~ stringr::str_remove_all(.x, "^\\[|\\]$")) %>%
  filter(`% IDY` > 87) %>%
  select(S1, E1) %>%
  arrange(S1)

merged <- iv(endpoint$S1, endpoint$E1 + 1) %>% iv_groups()
merged_df <- tibble(start = iv_start(merged), end = iv_end(merged) - 1)

# commentable section to bin coverage figure
binwidth <- nrow(cov) / 10000
tibble(s = seq(1, nrow(cov) - binwidth, binwidth),
       e = s + binwidth)

bin_starts <- seq(1L, nrow(cov) - binwidth, binwidth)

binned_cov <- tibble(
  Position = bin_starts + binwidth/2,
  Coverage = map2_dbl(bin_starts, bin_starts + binwidth, ~ median(cov$Coverage[.x:.y])))

# draw coverage plot
p1 <-
  binned_cov %>%
  #cov %>% # switch to this for no bin!
  ggplot(aes(x=Position, y=Coverage)) +
  geom_col(fill="#5D3A9B") +
  labs(
    x = NULL,
    y = "Fold10 Coverage"
  ) +
  theme_classic() +
  scale_x_continuous(labels = scales::comma, limits = c(1, max(cov$Position, merged_df$end))) +
  scale_y_continuous(trans='log10') + 
  theme(
    axis.title = element_text(size = 8),
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA)
  )


# draw contigs plot
p2 <- 
  ggplot(merged_df) +
  geom_segment(aes(x = start, xend = end, y = 0, yend = 0), colour = "#E66100", linewidth = 10) +
  annotation_custom(grid::textGrob("Contigs", x = grid::unit(0, "mm"), hjust = 0, gp = grid::gpar(fontsize = 8))) +
  labs(
    x = "Base position in genome",
    y = NULL
  ) +
  scale_x_continuous(labels = NULL, limits = c(1, max(cov$Position, merged_df$end))) +
  scale_y_continuous(labels = NULL) +
  theme_void() +
  theme(panel.background = element_rect(fill = "transparent",
                                            colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)
  )

# combine plots
p1 / p2 + 
  plot_layout(ncol = 1, heights = c(10, 1)) + 
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "transparent", color = NA)
    ))

# save!
ggsave(outname,
       width = 1500,
       height = 500,
       units = "px",
       dpi = 300,
       scale = 2,
       bg = "transparent")
