#!/usr/bin/env Rscript

# Load libraries
library(tidyverse)
library(ggplot2)
library(readr)
library(ggrepel)
library(plotly)
library(ggstar)


# This script reads a list of PCA result files (in `.evec` format), reads the prefix (without .evec).
# Generates customized PCA scatter plots for each, and saves them as PDF files.

# Assign unique colors, shapes, and transparency (alpha) to each category
colors = c(
  "Asia" = "#56B4E9", "Oceania" = "#FF69B4", "Turkey - Bezoar" = "#009E73",
  "Africa" = "#66A61E", "Other Capra" = "#FABEBE", "Europe" = "#E58A4D",
  "Iran - Bezoar" = "#666666", "Iran/Asia2" = "#9B1B30", "Russia" = "#7570B3",
  "0.1X" = "black", "0.15X" = "darkred", "0.5X" = '#0072B2',
  "0.25X" = "#FF7F00", "HQ" = "#Aaffc3", "0.75X" = "lightpink",
  "1X" = "turquoise", "0.3X" = "gray", "0.4X" = "#BCF60C",
  "2X" = "#000075", "4X" = "hotpink", "America" = '#301934'
)

shape = c(
  "Asia" = 15, "Oceania" = 15, "Turkey - Bezoar" = 23, "Africa" = 15,
  "Other Capra" = 14, "Europe" = 15, "Iran - Bezoar" = 23, "Iran/Asia2" = 15,
  "Russia" = 15, "Imputed" = 13, "HQ" = 28, "Haploid" = 11, "America" = 15
)

alpha = c(
  "Asia" = 0.41, "Oceania" = 0.41, "Turkey - Bezoar" = 0.41, "Africa" = 0.41,
  "Other Capra" = 0.41, "Europe" = 0.41, "Iran - Bezoar" = 0.41,
  "Iran/Asia2" = 0.41, "Russia" = 0.41, "Imputed" = 1, "HQ" = 1,
  "Haploid" = 1, "America" = 0.41
)

# Read the list of PCA file prefixes
files = read.table('PCA_list.txt')

# Loop through each PCA file

# Loop through each PCA file
for (x in files){
  for (id in x){

    # Read PCA .evec file (tab-delimited format)
    pca <- read.table(paste0(id, '.evec'), sep = "\t")
    # Load geographical, and downsample coverage information 
    info<-data.table("Vargoats_ancient_ID_origin_haploid_4col.txt", sep = "\t")
    pca_info<-left_join(pca, info, by = c('V1' = 'ID'),multiple = "any")

    # Recode population labels with defined order for plotting
    pca$style <- factor(pca$style, levels = c(
      'Europe', 'Russia', 'Africa', 'Asia', 'Oceania', 'America',
      'Iran/Asia2', 'Turkey - Bezoar', 'Iran - Bezoar',
      'Other Capra', 'Imputed', 'HQ', 'Haploid'))

    # Separate modern population samples
    Moderns <- pca[pca$style %in% c(
      'Asia', 'Oceania', 'Turkey - Bezoar', 'Iran - Bezoar', 'Africa',
      'Europe', 'Other Capra', 'Iran/Asia2', 'Russia'),]

    # Separate high-quality, imputed, and haploid samples
    Imp_HQ <- pca[pca$style %in% c("Imputed", "HQ", "Haploid"),]

    # Create the PCA plot
    plot_pca <- ggplot(pca, aes(x = V5, y = V6)) +
      geom_star(data = pca,
                aes(starshape = style, fill = col, alpha = style, colour = col),
                stat = "identity", size = 0.8, starstroke = 0.001) +
      labs(x = 'PC1', y = 'PC2', size = 3) +
      scale_fill_manual(values = colors) +
      scale_colour_manual(values = colors) +
      scale_starshape_manual(values = shape) +
      scale_alpha_manual(values = alpha) +
      geom_star(data = Imp_HQ,
                aes(x = V5, y = V6, starshape = style, fill = col),
                stat = "identity", size = 1.2, starstroke = 0.15, alpha = 1) +
      theme_minimal() +
      scale_x_reverse() +  
      scale_y_reverse() +
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = "black", size = 0.15),
        axis.text.y = element_text(size = 3),
        axis.text.x = element_text(size = 3),
        strip.text.x = element_text(size = 3),
        legend.position = 'right',
        legend.key.width = unit(0.12, 'cm'),
        legend.key.height = unit(0.025, 'cm'),
        legend.text = element_text(size = 3),
        legend.title = element_text("", size = 2),
        legend.spacing.y = unit(0.12, 'cm'),
        legend.spacing.x = unit(0.025, 'cm'),
        axis.ticks.y = element_line(size = 0.15),
        axis.ticks.x = element_line(size = 0.15),
        axis.title.y = element_text(size = 3),
        axis.title.x = element_text(size = 3)
      )

    # Export the PCA plot as a PDF
    pdf(paste0(id, '.pdf'), height = 2.4, width = 3.53)
    print(plot_pca)
    dev.off()
  }
}
