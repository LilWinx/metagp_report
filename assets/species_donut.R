list.of.packages <- c(
  "ggplot2",
  "tidyverse",
  "patchwork"
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
species <- read.csv("tpmAbundances.txt", header = FALSE, sep = "\t")
#outname <- "output_hbar.png"
outname <- paste(args[2],"_donut.png", sep = "")
dbPath <- args[3]
#dbPath <- "/Users/wfon4473/Documents/Bioinformatics/metagp_report/database"

# pre-processing for second hbar
fungi_db_path <- paste(dbPath,"/fungi.txt", sep = "")
fungi_db <- read.csv(fungi_db_path, header = 1, sep = "\t")
fungi_species <- fungi_db$X.Organism.Name

parasite_db_path <- paste(dbPath,"/parasites.txt", sep = "")
parasite_db <- read.csv(parasite_db_path, header = 1, sep = "\t")
parasite_species <- parasite_db$X.Organism.Name

kingdoms <- c("Bacteria", "Viruses")
selected <- species[species$Species %in% c(fungi_species, parasite_species) | species$Kingdom %in% kingdoms, ]
selected <- subset(selected, Phylum != "Chordata" & Genus != "Unknown")
s_total_TPM <- sum(selected$TPM) # total excluding Chordata & Unknowns
sorted <- selected[order(selected$TPM, decreasing = TRUE), ]
top_10 <- head(sorted, 10)
top_10 <- top_10[, c("Species", "TPM")]
top_10$Percentage <- (top_10$TPM / s_total_TPM) * 100
top10_total_tpm <- sum(top_10$Percentage)
other_tpm <- 100 - top10_total_tpm # calculate other
other_row <- data.frame(Species = "Other", TPM = 0, Percentage = other_tpm) # make new row
new_species <- rbind(top_10, other_row) # merge row into species data

colours <- c('#0079bf', 
             '#70b500', 
             '#ff9f1a', 
             '#eb5a46', 
             '#f2d600', 
             '#c377e0', 
             '#ff78cb', 
             '#00c2e0', 
             '#51e898',
             '#bc9b6a',
             '#c4c9cc')

species_freq <- table(new_species$Species)
sorted_species <- names(sort(species_freq, decreasing = TRUE))
sorted_species <- c(sorted_species[sorted_species != "Other"], "Other")
new_species$Species <- factor(new_species$Species, levels = sorted_species)

sorted <- mutate(new_species,
                 species = fct_infreq(new_species$Species, w = new_species$Percentage) |> fct_other(drop = "Other"))

donut <- ggplot(sorted, aes(x = 2, y = Percentage, fill=species)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.5, 2.5)) +
  scale_fill_manual(values = colours, name = "Species") +
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

ggsave(outname,
       width = 800,
       height = 500,
       units = "px",
       dpi = 300,
       scale = 3,
       bg = "transparent")
