list.of.packages <- c(
  "ggplot2",
  "tidyverse",
  "patchwork",
  "dplyr",
  "cowplot",
  "extrafont",
  "extrafontdb"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)
library(cowplot)
library(extrafont)

# load the fonts
loadfonts()

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
#wd <- "~/Documents/R_workdir"
setwd(dir = wd)
rna_species <- read.csv(args[2], header = 1)
#rna_species <- read.csv("~/Documents/R_workdir/CzRnaSoMuCSF_zscore.csv", header = 1)
dna_species <- read.csv(args[3], header = 1)
#dna_species <- read.csv("~/Documents/R_workdir/CzDnaSoMuCSF_zscore.csv", header = 1)
#outname <- "output_quaddonut.png"
outname <- paste(args[4],"_quaddonut.png", sep = "")
dbPath <- args[5]
#dbPath <- "/Users/wfon4473/Documents/Bioinformatics/metagp_report/database"
rna_species <- rna_species %>% filter(zscore > 1)

# set-up fungi database
fungi_db_path <- paste(dbPath,"/fungi.txt", sep = "")
fungi_db <- read.csv(fungi_db_path, header = 1, sep = "\t")
fungi_species <- fungi_db$X.Organism.Name

# set-up parasite database
parasite_db_path <- paste(dbPath,"/parasites.txt", sep = "")
parasite_db <- read.csv(parasite_db_path, header = 1, sep = "\t")
parasite_species <- parasite_db$X.Organism.Name

klist <- c(
  "Archaea",
  "Bacteria",
  "Eukaryota",
  "Viruses",
  "Other"
)

# pre-processing for first hbar
kingdom_process <- function(species) {
  filt_kingdom <- species[, c("superkingdom", "rpm_sample")] # get kingdom and rpm
  k_total_rpm <- sum(filt_kingdom$rpm_sample)
  filt_kingdom$superkingdom <- ifelse(filt_kingdom$superkingdom == "", "Other", filt_kingdom$superkingdom)
  sum_df <- aggregate(rpm_sample ~ superkingdom, data = filt_kingdom, FUN = sum)
  sum_df$percentage <- (sum_df$rpm_sample / k_total_rpm) * 100
  
  existing_kingdoms <- unique(sum_df$superkingdom)
  missing_kingdoms <- setdiff(klist, existing_kingdoms)
  if (length(missing_kingdoms) > 0) {
    missing_rows <- data.frame(superkingdom = missing_kingdoms, rpm_sample = 0, percentage = 0)
    sum_df <- rbind(sum_df, missing_rows)
  }
  
  kingdom_freq <- table(sum_df$superkingdom)
  sorted_kingdom <- names(sort(kingdom_freq, decreasing = TRUE))
  sorted_kingdom <- c(sorted_kingdom[sorted_kingdom != "Other"], "Other")
  sum_df$superkingdom <- factor(sum_df$superkingdom, levels = sorted_kingdom)
  return(sum_df)
}



# pre-processing for second hbar
top10species <- function(species) {
  kingdoms <- c("Bacteria", "Viruses")
  selected <- species[species$species %in% c(fungi_species, parasite_species) | species$superkingdom %in% kingdoms, ]
  selected <- subset(selected, phylum != "Chordata" & genus != "Unknown")
  s_total_rpm <- sum(species$rpm_sample) # total excluding Chordata & Unknowns
  sorted <- subset(selected, !is.na(species) & species != "")
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

  species_freq <- table(new_species$species)
  sorted_species <- names(sort(species_freq, decreasing = TRUE))
  sorted_species <- c(sorted_species[sorted_species != "Other"], "Other")
  new_species$species <- factor(new_species$species, levels = sorted_species)

  sorted <- mutate(new_species,
                  species = fct_infreq(new_species$species, w = new_species$percentage) |> fct_other(drop = "Other"))
                  
  return(sorted)
}

# COLOURS
kcolours <- c(
  Archaea = '#ef476f',
  Bacteria = '#ffd166',
  Eukaryota = '#06d6a0',
  Viruses = '#118ab2',
  Other = '#c4c9cc'
)



colours <- c('#0079bf', 
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

generate_colors <- function(n) {
  colours <- c('#0079bf', 
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
  
  if (n <= 0) {
    return(c(Other = '#c4c9cc'))
  } else if (n <= length(colours)) {
    return(c(colours[1:n], Other = '#c4c9cc'))
  } else {
    dynamic_colors <- rainbow(n - 1)
    return(c(dynamic_colors, Other = '#c4c9cc'))
  }
}

# Draw the kingdom plot
rna_kingdom_df <- kingdom_process(rna_species)
dna_kingdom_df <- kingdom_process(dna_species)

kingdom_donut <- function(kingdom_df, na_type) { 
  donut <- ggplot(kingdom_df, aes(x = 2, y = percentage, fill=superkingdom)) +
    geom_col() +
    coord_polar(theta = "y") +
    xlim(c(0.5, 2.5)) +
    scale_fill_manual(values = kcolours, name = "Kingdom", limits = klist, drop = FALSE) +
    theme_void() +
    guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7)) + 
    theme(
      legend.position="bottom",
      legend.direction='vertical',
      legend.text = element_text(size = 7, family = "DejaVu Sans"), 
      legend.title = element_text(size = 9, family = "DejaVu Sans"),
      panel.background = element_rect(fill = "transparent",
                                      colour = NA_character_),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent",
                                     colour = NA_character_),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(b = 2.5, unit = "pt")
    )
  return(donut)  
}

rna_kingdom_donut <- kingdom_donut(rna_kingdom_df)
dna_kingdom_donut <- kingdom_donut(dna_kingdom_df)

# Draw the species plot
rna_species_df <- top10species(rna_species)
rna_n <- nrow(rna_species_df) - 1
dna_species_df <- top10species(dna_species)
dna_n <- nrow(dna_species_df) - 1

rna_colours <- generate_colors(rna_n)
dna_colours <- generate_colors(dna_n)

species_donut <- function(species_df, na_type, na_colours) { donut <- ggplot(species_df, aes(x = 2, y = percentage, fill=species)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.5, 2.5)) +
  scale_fill_manual(values = setNames(na_colours, levels(species_df$species)), name = "Species") +
  theme_void() +
  guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7)) +
  theme(
    legend.position="right", 
      legend.text = element_text(size = 7, family = "DejaVu Sans"), 
      legend.title = element_text(size = 9, family="DejaVu Sans"),
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
    ) +
    annotate("text", x = 0.5, y = 0.5, label = na_type, size = 7, family = "DejaVu Sans")
  return(donut)  
  }

rna_species_donut <- species_donut(rna_species_df, "RNA", rna_colours)
dna_species_donut <- species_donut(dna_species_df, "DNA", dna_colours)

# combine donuts
kingdom <- (rna_kingdom_donut + theme(plot.margin = unit(c(0,0,70,0), "pt"))) / (dna_kingdom_donut + theme(plot.margin = unit(c(0,0,0,0), "pt")))+
  plot_layout(ncol = 1, heights = c(0.3,0.7), guides = "collect") + 
  plot_annotation(theme = theme(
        plot.background = element_rect(fill = "transparent", color = NA)
      )) & theme(legend.position = 'bottom')
species <- (rna_species_donut + theme(plot.margin = unit(c(0,0,0,0), "pt"))) / (dna_species_donut + theme(plot.margin = unit(c(0,0,0,0), "pt")))+
  plot_layout(ncol = 1, heights = c(3,3)) +
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "transparent", color = NA)
  )) & theme(legend.position = 'right')

p <- plot_grid(kingdom, NULL, species, ncol = 3, nrow = 1, rel_widths = c(1, -0.5, 4))

ggsave(outname,
       width = 1450,
       height = 1400,
       units = "px",
       dpi = 300,
       scale = 1.2,
       bg = "transparent")

