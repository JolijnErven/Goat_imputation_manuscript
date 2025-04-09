#!/usr/bin/env Rscript

# Load libraries
require(ggplot2)
library(tidyverse)
library(rnaturalearth)
library(sf)
library(raster)
library(ggrepel)
library(wesanderson)
require(ggnewscale)

# Script to plot outgroup f3 statistics -- created by Kevin Daly
# This script generates geographic visualizations of outgroup f3 statistics.
# It includes:
#   - Ancient-to-Modern comparisons (modern goat breeds vs ancient genomes)
#   - Ancient-to-Ancient comparisons
# Data is pulled from local files and supplementary tables hosted via Google Sheets.
# Script created by Kevin Daly


# Load geographical map
land <- ne_load(scale = 50, type = "land", category = "physical", destdir = "/home/kdaly/Downloads/ne_50m_land/", returnclass = "sf")
countries <- ne_countries( returnclass = "sf", scale = 50)

# Reads ancient to ancient F3 out group statistics
f3.ancient<-read.table("vargoats_f3-out_march2025_combined_ancients-ancients.out",col.names = c("A","B","Out","f3","err","Z"))
# Reads ancient to modern F3 out group statistics
f3.modern<-read.table("vargoats_f3-out_march2025_combined_ancients-moderns.out",col.names = c("A","B","Out","f3","err","Z"))

# Get ancient metadatafrom Table S3
anc.coord<-data.frame(read_sheet("https://docs.google.com/spreadsheets/d/1og33Y8H9PMYM2DTgDoM9vN63fRROa7v1E21twzl3UJc/edit?gid=1991801573#gid=1991801573",sheet = 2))
anc.coord <- data.frame(as.matrix(anc.coord))
anc.coord$Sample<-unlist(anc.coord$Sample)
# Get modern metadata from Table S1
vargoats.coord<-data.frame(read_sheet("https://docs.google.com/spreadsheets/d/1nF7_hAWWKhrW7K2xVYj1ci2HUoQdEE56qegfHcEE6xM/edit?gid=1922871814#gid=1922871814"))
vargoats.coord <- data.frame(as.matrix(vargoats.coord))
vargoats.coord$combined_new_breed_code_priority_for_split_off<-unlist(vargoats.coord$combined_new_breed_code_priority_for_split_off)
vargoats.coord.sub<-subset(vargoats.coord,select=c("combined_new_breed_code_priority_for_split_off","Latitude.disctric","Longitude.country","Latitude.country","Longitude.district","Species"))

# Merge Modern outgroup f3 statistics with metadata
f3.modern.merge<-left_join(f3.modern, vargoats.coord.sub, by = c('A' = 'combined_new_breed_code_priority_for_split_off'),multiple = "any")

# Transform to numeric
f3.modern.merge$Longitude.country<-as.numeric(as.character(f3.modern.merge$Longitude.country))
f3.modern.merge$Longitude.district<-as.numeric(as.character(f3.modern.merge$Longitude.district))
f3.modern.merge$Latitude.country<-as.numeric(as.character(f3.modern.merge$Latitude.country))
f3.modern.merge$Latitude.disctric<-as.numeric(as.character(f3.modern.merge$Latitude.disctric))
anc.coord$long<-as.numeric(as.character(anc.coord$long))
anc.coord$lat<-as.numeric(as.character(anc.coord$lat))

# Split by species: Domestic vs Wild Capra
f3.modern.merge.CH<-f3.modern.merge[f3.modern.merge$Species=="Capra hircus",]
f3.modern.merge.Capra<-f3.modern.merge[f3.modern.merge$Species!="Capra hircus",]
f3.modern.merge.Capra$Species<-unlist(f3.modern.merge.Capra$Species)

# Create geographical map
blankmap <- ggplot() + geom_sf(data = land, color = "grey50", fill = "aliceblue", linewidth=1) + 
  geom_sf(data = countries, color = "white", fill = "white", linewidth = 1.0) + 
  theme(panel.grid.major = element_line(color = "aliceblue", linewidth = 0), 
        panel.background = element_rect(fill = "aliceblue")) + 
  scale_shape_manual(values = c(22,21)) + theme(legend.position = "right") + 
  theme(axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
  theme(axis.title.y = element_text(size=14)) +theme(axis.title.x = element_text(size=14)) +
  xlim(-80,175) + ylim(c(-55,80))

# Colour pallete
pal<-wes_palette("Zissou1", n = 100, type="continuous")

# Code domestic breeds as domestic
f3.modern.merge.CH$species.coding<-"Domestic"

# Generate Modern vs Ancient f3 Maps; Loop through every ancient sample
for (i in unique(f3.modern.merge.CH$B)) {
a<-blankmap + new_scale_fill() + 
  geom_jitter(data = f3.modern.merge.CH[f3.modern.merge.CH$B==i,], alpha=0.95, shape=21, 
              size=4, aes(y=Latitude.disctric, x=Longitude.district, fill=f3),width=2,height = 2)  +
  scale_fill_gradientn(colors = pal)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Shared drift between modern breed grouping and ",i)) +
  xlab("") + ylab("")
svg(paste0("shared_drift_",i,"_vargoats.svg"),width=10, height = 9)
print(a)
dev.off()
}

# Plotting ancient to ancient f3s
# Merge ancient f3 statistics with information
f3.ancient.merge<-left_join(f3.ancient, anc.coord, by = c('A' = 'name'),multiple = "any")

# Create geographical map for ancients
blankmap.anc <- ggplot() + geom_sf(data = land, color = "grey50", fill = "aliceblue", linewidth=1) + 
  geom_sf(data = countries, color = "white", fill = "white", linewidth = 1.0) + 
  theme(panel.grid.major = element_line(color = "aliceblue", linewidth = 0), 
        panel.background = element_rect(fill = "aliceblue")) + 
  scale_shape_manual(values = c(22,21)) + theme(legend.position = "right") + 
  theme(axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
  theme(axis.title.y = element_text(size=14)) +theme(axis.title.x = element_text(size=14)) +
  xlim(-5,72) + ylim(c(27.5,55.5))

# For every ancient create a plot
for (i in unique(f3.ancient.merge$B)) {
  a<-blankmap.anc + new_scale_fill() + 
  geom_jitter(data = f3.ancient.merge[f3.ancient.merge$B==i,], alpha=0.95, 
              size=4, aes(y=lat, x=long, fill=f3,shape=period.plotting.expanded),width=0.5,height = 0.5)  +
    scale_shape_manual(values=c(21,22,24,23,25)) +
  scale_fill_gradientn(colors = pal)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())  +
    ggtitle(paste0("Shared drift between ancient genome and ",i)) +
    xlab("") + ylab("")
  svg(paste0("shared_drift_",i,"_ancients.svg"),width=10, height = 9)
  print(a)
  dev.off()
}



