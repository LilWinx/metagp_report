library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)

#wd <- args[1]
wd <- "/Users/wfon4473/Documents/R_workdir"
setwd(dir = wd)
species <- read.csv("tpmAbundances.txt", header = FALSE, sep = "\t")
outname <- paste("test","_hbar.png", sep = "")
#outname <- "test.png"

# pre-processing for first hbar
colnames(species) <- c("TPM","Kingdom","Phylum","Class","Order","Family","Genus","Species")
filt_kingdom <- species[, c("Kingdom", "TPM")] # get kingdom and TPM
k_total_TPM <- sum(filt_kingdom$TPM)
filt_kingdom$Kingdom <- ifelse(filt_kingdom$Kingdom == "Unknown", "Other", filt_kingdom$Kingdom)
sum_df <- aggregate(TPM ~ Kingdom, data = filt_kingdom, FUN = sum)
sum_df$Percentage <- (sum_df$TPM / k_total_TPM) * 100

kingdom_freq <- table(sum_df$Kingdom)
sorted_kingdom <- names(sort(kingdom_freq, decreasing = TRUE))
sorted_kingdom <- c(sorted_kingdom[sorted_kingdom != "Other"], "Other")
sum_df$Kingdom <- factor(sum_df$Kingdom, levels = sorted_kingdom)

# pre-processing for second hbar
fungi_phyla <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Microsporidia", "Glomeromycota", "Zygomycota")
parasite_phyla <- c("Apicomplexa", "Ciliophora", "Bacillariophyta", "Cercozoa", "Euglenozoa", "Heterolobosea", "Parabasalia", "Fornicata", "Evosea", "Streptophyta")
kingdoms <- c("Bacteria", "Viruses")
selected <- species[species$Phylum %in% c(fungi_phyla, parasite_phyla) | species$Kingdom %in% kingdoms, ]
selected <- subset(selected, Phylum != "Chordata")
s_total_TPM <- sum(selected$TPM) # total excluding Chordata.
sorted <- selected[order(selected$TPM, decreasing = TRUE), ]
top_10 <- head(sorted, 10)
top_10 <- top_10[, c("Species", "TPM")]
top_10$Percentage <- (top_10$TPM / s_total_TPM) * 100
top10_total_tpm <- sum(top_10$Percentage)
other_tpm <- 100 - top10_total_tpm # calculate other
other_row <- data.frame(Species = "Other", TPM = 0, Percentage = other_tpm) # make new row
new_species <- rbind(top_10, other_row) # merge row into species data

kcolours <- c(
  '#ef476f',
  '#ffd166',
  '#06d6a0',
  '#118ab2',
  '#c4c9cc'
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
             '#bc9b6a',
             '#c4c9cc')

species_freq <- table(new_species$Species)
sorted_species <- names(sort(species_freq, decreasing = TRUE))
sorted_species <- c(sorted_species[sorted_species != "Other"], "Other")
new_species$Species <- factor(new_species$Species, levels = sorted_species)

sorted <- mutate(new_species,
                 species = fct_infreq(new_species$Species, w = new_species$Percentage) |> fct_other(drop = "Other"))

hbar1 <- ggplot(sum_df, aes(x = 2, y = Percentage, fill=Kingdom)) +
  geom_col(width = 0.1) +
  coord_flip() + 
  scale_y_reverse() +
  scale_fill_manual(values = kcolours, name = "Kingdom") +
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
        legend.key = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(b = 2.5, unit = "pt")
  )

hbar2 <- ggplot(sorted, aes(x = 2, y = Percentage, fill=Species)) +
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
