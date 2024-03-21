list.of.packages <- c(
  "ggplot2",
  "tidyverse",
  "patchwork",
  "dplyr"
  )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
#wd <- "~/Documents/R_workdir"
setwd(dir = wd)
files <- list.files(wd, pattern = "_zscore\\.csv$", full.names = TRUE)
if (length(files) == 1) {
    # Assign the file path to a variable
    species <- read.csv(files[1], header = 1)
} else {
    print("Error: No file or multiple files found with '_zscore.csv' ending.")
}
#outname <- "output_hbar.png"
outname <- paste(args[2],"_hbar.png", sep = "")
dbPath <- args[3]
#dbPath <- "/Users/wfon4473/Documents/Bioinformatics/metagp_report/database"
species <- species %>% filter(zscore > 1)

# pre-processing for first hbar
filt_kingdom <- species[, c("superkingdom", "rpm_sample")] # get kingdom and rpm
k_total_rpm <- sum(filt_kingdom$rpm_sample)
filt_kingdom$superkingdom <- ifelse(filt_kingdom$superkingdom == "", "Other", filt_kingdom$superkingdom)
sum_df <- aggregate(rpm_sample ~ superkingdom, data = filt_kingdom, FUN = sum)
sum_df$percentage <- (sum_df$rpm_sample / k_total_rpm) * 100

kingdom_freq <- table(sum_df$superkingdom)
sorted_kingdom <- names(sort(kingdom_freq, decreasing = TRUE))
sorted_kingdom <- c(sorted_kingdom[sorted_kingdom != "Other"], "Other")
sum_df$superkingdom <- factor(sum_df$superkingdom, levels = sorted_kingdom)

# pre-processing for second hbar
fungi_db_path <- paste(dbPath,"/fungi.txt", sep = "")
fungi_db <- read.csv(fungi_db_path, header = 1, sep = "\t")
fungi_species <- fungi_db$X.Organism.Name

parasite_db_path <- paste(dbPath,"/parasites.txt", sep = "")
parasite_db <- read.csv(parasite_db_path, header = 1, sep = "\t")
parasite_species <- parasite_db$X.Organism.Name



kingdoms <- c("Bacteria", "Viruses")
selected <- species[species$species %in% c(fungi_species, parasite_species) | species$superkingdom %in% kingdoms, ]
selected <- subset(selected, phylum != "Chordata" & genus != "Unknown")
s_total_rpm <- sum(species$rpm_sample) # total excluding Chordata & Unknowns
sorted <- selected[order(selected$rpm_sample, decreasing = TRUE), ]
sorted <- subset(sorted, !is.na(species) & species != "")
sorted <- sorted %>%
  group_by(species) %>%
  summarize(rpm_sample = sum(rpm_sample), zscore = max(zscore)) %>%
  ungroup()
sorted <- sorted[order(sorted$zscore, decreasing = TRUE), ]
top_10 <- sorted %>%
  filter(zscore > 50) %>%
  head(10)
top_10 <- top_10[, c("species", "rpm_sample")]

top_10$percentage <- (top_10$rpm_sample / s_total_rpm) * 100
top10_total_rpm <- sum(top_10$percentage)
other_rpm <- 100 - top10_total_rpm # calculate other
other_row <- data.frame(species = "Other", rpm_sample = 0, percentage = other_rpm) # make new row
top_10 <- top_10[, c("species", "rpm_sample", "percentage")]
new_species <- rbind(top_10, other_row) # merge row into species data


kcolours <- c(
  '#ef476f',
  '#ffd166',
  '#06d6a0',
  '#118ab2',
  '#c4c9cc'
)

generate_colors <- function(n) {
  if (n <= 0) {
    return(c(Other = '#c4c9cc'))
  } else if (n == 1) {
    return(c('#0079bf', '#c4c9cc'))  # Return a default color when only one color is needed
  } else {
    colours <- c(
      '#0079bf', 
      '#70b500', 
      '#ff9f1a', 
      '#eb5a46', 
      '#f2d600', 
      '#c377e0', 
      '#ff78cb', 
      '#00c2e0', 
      '#51e898',
      '#bc9b6a'
    )
    if (n <= length(colours)) {
      return(c(colours[1:n], '#c4c9cc'))
    } else {
      dynamic_colors <- rainbow(n - 1)
      return(c(dynamic_colors, '#c4c9cc'))
    }
  }
}


species_freq <- table(new_species$species)
sorted_species <- names(sort(species_freq, decreasing = TRUE))
sorted_species <- c(sorted_species[sorted_species != "Other"], "Other")
new_species$species <- factor(new_species$species, levels = sorted_species)

sorted <- mutate(new_species,
                 species = fct_infreq(new_species$species, w = new_species$percentage) |> fct_other(drop = "Other"))

sorted_n <- nrow(sorted) - 1
hbar2_colours <- generate_colors(sorted_n) 
#hbar2_colours <- c("#0079bf","#c4c9cc")

hbar1 <- ggplot(sum_df, aes(x = 2, y = percentage, fill=superkingdom)) +
  geom_col(width = 0.1) +
  coord_flip() + 
  scale_y_reverse() +
  scale_fill_manual(values = kcolours, name = "Kingdom") +
  theme_void() +
  guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="bottom", 
        legend.text = element_text(size = 7, family = "Arial"), 
        legend.title = element_text(size = 9, family="Arial"),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(b = 2.5, unit = "pt")
  )

hbar2 <- ggplot(sorted, aes(x = 2, y = percentage, fill=species)) +
                geom_col(width = 0.1) +
                coord_flip() + 
                scale_y_reverse() +
                scale_fill_manual(values = hbar2_colours, name = "Species") +
                theme_void() +
                guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7)) +
                theme(legend.position="bottom", 
                      legend.text = element_text(size = 7, family="Arial"), 
                      legend.title = element_text(size = 9, family="Arial"),
                      panel.background = element_rect(fill = "transparent",
                                                      colour = NA_character_), # necessary to avoid drawing panel outline
                      panel.grid.major = element_blank(), # get rid of major grid
                      panel.grid.minor = element_blank(), # get rid of minor grid
                      plot.background = element_rect(fill = "transparent",
                                                     colour = NA_character_), # necessary to avoid drawing plot outline
                      legend.background = element_rect(fill = "transparent", color = NA),
                      legend.box.background = element_rect(fill = "transparent", color = NA),
                      legend.key = element_rect(fill = "transparent", color = NA),
                      plot.margin = margin(t = 2.5, unit = "pt")
                )
# combine hbars
p <- hbar1 / hbar2 + 
        plot_layout(ncol = 1, heights = c(1, 1)) + 
        plot_annotation(theme = theme(
          plot.background = element_rect(fill = "transparent", color = NA)
        )) 

ggsave(outname,
       width = 1500,
       height = 300,
       units = "px",
       dpi = 300,
       scale = 1.7,
       bg = "transparent")
