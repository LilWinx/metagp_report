library(ggplot2)
library(tidyverse)
library(dplyr)

#setwd(dir = "/Users/wfon4473/Documents/R_workdir")
species <- read_tsv("species.txt")
outname <- "horizontalbar.png"

# pre-processing
filt_species <- species[, c("species", "total_RA")] # get only the two columns needed
total_abundance <- sum(filt_species$total_RA) # calculate total_RA from Kraken
other_abundance <- 100 - total_abundance # calculate other
other_row <- data.frame(species = "Other", total_RA = other_abundance) # make new row
new_species <- rbind(filt_species, other_row) # merge row into species data
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

species_freq <- table(new_species$species)
sorted_species <- names(sort(species_freq, decreasing = TRUE))
sorted_species <- c(sorted_species[sorted_species != "Other"], "Other")
new_species$species <- factor(new_species$species, levels = sorted_species)

sorted <- mutate(new_species,
                 species = fct_infreq(new_species$species, w = new_species$total_RA) |> fct_other(drop = "Other"))

ggplot(sorted, aes(x = 2, y = total_RA, fill=species)) +
  geom_col(width = 0.1) +
  coord_flip() + 
  scale_y_reverse() +
  scale_fill_manual(values = colours, name = "Species") +
  theme_void() +
  guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="bottom", 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9),
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

ggsave(outname,
       width = 1500,
       height = 150,
       units = "px",
       dpi = 300,
       scale = 1.7,
       bg = "transparent")
