#!/usr/bin/env Rscript

# Load libraries
library(tidyverse)
library(ggthemes)
library(googlesheets4)  
library("rnaturalearth")
library("rnaturalearthdata")
library(sf)
library(ggrepel)
library(RColorBrewer)

# Script to plot IBD -- creates Figure 6

# Load necessary files for geographical map
land <- ne_download(
  scale = 50,
  type = "land",
  category = "physical",
  returnclass = "sf")

rivers <- ne_download(scale=50, type='rivers_lake_centerlines',category = "physical", returnclass = "sf")
lakes <- ne_download(scale=50, type='lakes',category = "physical", returnclass = "sf")
world1 <- ne_countries( returnclass = "sf",scale = 50)


# Load dataframe from Table S15 Normalized relatedness -- creates Figure 6A
IBD_df <-data.frame(read_sheet("https://docs.google.com/spreadsheets/d/1og33Y8H9PMYM2DTgDoM9vN63fRROa7v1E21twzl3UJc/edit?gid=1991801573#gid=1991801573",sheet = 15))
# Transform longitute, latitude and normalized relatedness as numeric
data_ibd_nr <- data_ibd %>% mutate_at(c('Long', 'Lat','normalized.relatedness'), as.numeric)

# Get colours
colors <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(10))

# Plot
plot_normIBD <- ggplot() +
  geom_sf(data=land, color = "darkgray", fill = "aliceblue", linewidth=0.15) + 
  geom_sf(data = world1, color = "darkgray", fill = "#F5F5F5", size = 0.053, alpha=0.6) +
  geom_point(data= data_ibd_nr,  mapping = aes(x=Long, y=Lat,fill=normalized.relatedness), shape=21,stroke =0.25, size =2.5) +
  geom_text_repel(data=data_ibd_nr,aes(label=as.factor(settlement)), colour='black') +
  scale_shape_manual(values = shapes ) +
  coord_sf(ylim=c(25,50), xlim=c(0,60)) +
  xlab("") + ylab("") +
  scale_fill_gradientn(colours=colors) +
  theme_classic(base_size = 3) +
  theme(panel.background = element_rect(fill='aliceblue'),
        panel.grid.major = element_line(color = "aliceblue", size = 0),
        axis.line = element_line(colour = 'black'),
        panel.border = element_rect(colour='black', fill=NA, size=0.5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'left'
  )

pdf("Final_Geographic_map_IBD_norm_relatedness_goats_Figure6A.pdf",width=3,height=3)
plot_normIBD
dev.off()

# Between settlement IBD cumulative cM -- creates Figure 6B
# Create dataframe with cumulative sum of cM between settlement based on Table S16 IBD between pairs
IBD_between <- data.frame(
  ID1 = c("Monjukli Depe", "Ganj Dareh", "Ganj Dareh", "Sang-e Chakhmaq", "Sang-e Chakhmaq", "Tepe Abdul Hosein"),
  Segments = c(3, 1, 1, 3, 1, 1),
  cumulative_cM = c(33.83445, 9.84972, 6.89988, 33.83445, 9.84972, 6.89988),
  pair = c(1, 2, 3, 1, 2, 3),
  Lat = c(36.8485, 34.27218, 34.27218, 36.504158, 36.504158, 34.066),
  Long = c(60.4181, 47.476423, 47.476423, 55.000786, 55.000786, 48.134)
)

# Deduplicate settlements by coordinates to avoid label overlap
unique_labels <- IBD_between %>%
  distinct(Lat, Long, .keep_all = TRUE)


plot_betweenIBD <- ggplot() +
  geom_sf(data=land, color = "black", fill = "aliceblue", linewidth=0.005) + 
  geom_sf(data = world1, color = "black", fill = "#F5F5F5", size = 0.0053, alpha=0.6) +
  geom_line(data= IBD_between,  mapping = aes(x=Long, y=Lat, color=cumulative_cM,group=pair), size =1) +
  geom_point(data= IBD_between,  mapping = aes(x=Long, y=Lat),stroke =0.05, size =1.5) +
  geom_text_repel(data = unique_labels, aes(x = Long, y = Lat, label = ID1), colour = 'black', size = 2) +
  scale_colour_gradientn(colours=colors, limits = c(0, 35),oob = scales::squish) +
  coord_sf(xlim=c(45,62), ylim=c(32,38)) +
  xlab("") + ylab("") +
  theme_classic(base_size = 5) +
  theme(panel.background = element_rect(fill='aliceblue'),
        panel.grid.major = element_line(color = "aliceblue", size = 0),
        axis.line = element_line(colour = 'black'),
        panel.border = element_rect(colour='black', fill=NA, size=0.25),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
  )

pdf("Neolithic_between_settleements_GP99_below_20_missing_glimpse_MAF_0.05_TV.LOD3_3cM_Figure6B.pdf",width=3,height=3)
plot_betweenIBD
dev.off()
