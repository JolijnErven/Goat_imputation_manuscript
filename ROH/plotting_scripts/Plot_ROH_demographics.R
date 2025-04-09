#!/usr/bin/env Rscript

# Load librarieslibrary(optparse)
library(tidyverse)
library(ggplot2)
library(forcats)
library(ggthemes)
library(grid)

# This R script generates demographic ROH (Runs of Homozygosity) bar plots 
# for figures used in the paper, this will be noted down at the start of a 
# new section


# Load dataframe from output Plink_ROH_pipeline_SNP_exploration.sh script -- Creates figure S28
ROH_file <-read.table("Vargoats_MAF_0.05_ancients_0.5x_HQ_similar_pos_no_low_GP_no_low_conc_no_missing_anc_1het_SNPs_comparison.txt", sep="\t",header=TRUE)

# Change downsampled_sites to easier to read variables and then sort
ROH_file$Downsampled_sites <- gsub("500000", "500k sites", ROH_file$Downsampled_sites)
ROH_file$Downsampled_sites <- gsub("750000", "750k sites", ROH_file$Downsampled_sites)
ROH_file$Downsampled_sites <- gsub("1000000", "1M sites", ROH_file$Downsampled_sites)
ROH_file$Downsampled_sites <- factor(ROH_file$Downsampled_sites,levels =c('Full',"1M sites",'750k sites',"500k sites"))


# Convert ROH length from kb to Mb
ROH_file$ROH_SUM <- as.numeric(ROH_file$ROH_SUM)/1000

# Set colours for ROH bins
colors <- c("#800000","#C21807","#FF817E","#D1EAF0","#0077B6","#012A4A")

# Plot
plot_roh <- ROH_file %>%
  mutate(Bins=fct_relevel(Bins,'>16.0','8.0-16.0','4.0-8.0','2.0-4.0','1.0-2.0','0.5-1.0')) %>% 
  ggplot(aes(fill=Bins, x=Origin, y=ROH_SUM)) +
  geom_col(colour='black', size = 0.3,width = 0.55) +
  scale_fill_manual(values = colors, name="ROH bins") +
  labs(x= '', y = "ROH length (Mb)", size = 4) +
  facet_grid(interaction(Sites,Dataset) ~ Downsampled_sites) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,350)) +
  scale_x_discrete(expand = expand_scale(add = 0.3)) +
  theme_minimal() +
  theme(
    panel.margin = grid::unit(0.5, "lines"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major.y = element_line( size=.03, color="lightgray"),
    panel.grid.minor.y = element_line( size=.03, color="lightgray"),
    axis.text.y = element_text(size=4),
    axis.text.x = element_text(angle = 90,size=4),
    strip.text.x = element_text(size = 4),
    legend.key.width= unit(0.35, 'cm'),
    legend.key.height= unit(0.35, 'cm'),
    legend.text = element_text(size=5),
    legend.title = element_text(size=5),
    axis.title.y = element_text(size = 5)
  ) 

# Plot to PDF
pdf("Vargoats_MAF_0.05_ancients_0.5x_HQ_similar_pos_no_low_GP_no_low_conc_no_missing_anc_1het_SNPs_comparison.pdf",width=9,height=5)
print(plot_roh)
dev.off()


# Load dataframe from output ROH_calculation_individual_no_missingness.sh script -- Creates figure S29.
ROH_file <-read.table("ROH_demo_plots_MAF_0.05_ancients_similar_pos_no_low_GP_no_low_conc_no_missing_500kb_100kb_200SNPs_1het_sum_transpose.txt", sep="\t",header=TRUE)

ROH_file$value <- as.numeric(ROH_file$value)/1000

Imp_HQ <- ROH_file[ROH_file$Origin %in% c("0-5X","HQ","0-25X","0-75X","1X"),]

colors <- c("#800000","#C21807","#FF817E","#D1EAF0","#0077B6","#012A4A")

plot_roh <- ROH_file %>%
  mutate(Bins=fct_relevel(Bins,'>16.0 Mb','8.0-16.0 Mb','4.0-8.0 Mb','2.0-4.0 Mb','1.0-2.0 Mb','0.5-1.0 Mb')) %>% 
  ggplot(aes(fill=Bins, x=Origin, y=ROH_SUM)) +
  geom_col(colour='black', size = 0.3,width = 0.55) +
  scale_fill_manual(values = colors, name="ROH bins") +
  labs(x= '', y = "ROH length (Mb)", size = 4) +
  facet_grid(~factor(ROH_file$Sample),scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0),limits = c(0,400)) +
  scale_x_discrete(expand = expand_scale(add = 0.6)) +
  theme_minimal() +
  theme(
    panel.margin = grid::unit(0.5, "lines"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major.y = element_line( size=.03, color="lightgray"),
    panel.grid.minor.y = element_line( size=.03, color="lightgray"),
    axis.text.y = element_text(size=4),
    axis.text.x = element_text(angle = 90,size=4),
    strip.text.x = element_text(size = 4),
    legend.key.width= unit(0.35, 'cm'),
    legend.key.height= unit(0.35, 'cm'),
    legend.text = element_text(size=5),
    legend.title = element_text(size=5),
    axis.title.y = element_text(size = 5)
  ) 

pdf("Final_ROH_demo_plots_MAF_0.05_ancients_similar_pos_no_low_GP_no_low_conc_no_missing_500kb_100kb_200SNPs_1het_sum.pdf",width=4,height=2)
print(plot_roh)
dev.off()


# Load dataframe from output BCFtools_ROH_no_missingness_imputed_test_samples.sh script -- Creates figure S30.
ROH_file <-read.table("ROH_demo_plotssort_vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_MAF_0.05_ancients_0.5x_10Q_200SNPs_sum_transpose.txt", sep="\t",header=TRUE)

# Note different from Plink, BCFtools is in bp so divide by 1000000 to get to Mb
ROH_file$ROH_SUM <- as.numeric(ROH_file$ROH_SUM)/1000000

colors <- c("#800000","#C21807","#FF817E","#D1EAF0","#0077B6","#012A4A")

plot_roh <- ROH_file %>%
  mutate(Bins=fct_relevel(Bins,'>16.0 Mb','8.0-16.0 Mb','4.0-8.0 Mb','2.0-4.0 Mb','1.0-2.0 Mb','0.5-1.0 Mb')) %>% 
  ggplot(aes(fill=Bins, x=Origin, y=ROH_SUM)) +
  geom_col(colour='black', size = 0.3,width = 0.55) +
  scale_fill_manual(values = colors, name="ROH bins") +
  labs(x= '', y = "ROH length (Mb)", size = 4) +
  scale_y_continuous(expand = c(0, 0.1),limits = c(0,400)) +
  scale_x_discrete(expand = expand_scale(add = 0.3)) +
  theme_minimal() +
  theme(
    panel.margin = grid::unit(0.5, "lines"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major.y = element_line( size=.03, color="lightgray"),
    panel.grid.minor.y = element_line( size=.03, color="lightgray"),
    axis.text.y = element_text(size=4),
    axis.text.x = element_text(angle = 90,size=4),
    strip.text.x = element_text(size = 4),
    legend.key.width= unit(0.35, 'cm'),
    legend.key.height= unit(0.35, 'cm'),
    legend.text = element_text(size=5),
    legend.title = element_text(size=5),
    axis.title.y = element_text(size = 5)
  ) 

pdf("ROH_demo_plotssort_vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_MAF_0.05_ancients_0.5x_10Q_200SNPs_sum_bcftools.pdf",width=4,height=2)
print(plot_roh)
dev.off()

# Load dataframes from output ROH_ancient_imputed_paleogenomes_PLINK_BCFtools.sh script -- Creates figure S32 and OSF figure.
ROH_file1 <- read.table("ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness_500kb_100kb_1het_200SNPs_sum_transpose_info.txt",  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ROH_file2 <- read.table("Combined_imputed_ancients_MAF_0.05_TV_1.2M_500kb_100kb_200SNPs_1het_sum_transpose_info.txt",  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ROH_file3 <- read.table("Combined_imputed_ancients_MAF_0.05_TV_1.6M_500kb_100kb_200SNPs_1het_sum_transpose_info.txt",  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ROH_file4 <- read.table("Combined_imputed_ancients_MAF_0.05_TV_1.6M_200SNPs_10Q_sum_transpose_info.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ROH_file5 <- read.table("ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q_sum_transpose_info.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Transform to MB
ROH_file1$ROH_SUM <- as.numeric(ROH_file1$ROH_SUM)/1000
ROH_file2$ROH_SUM <- as.numeric(ROH_file2$ROH_SUM)/1000
ROH_file3$ROH_SUM <- as.numeric(ROH_file3$ROH_SUM)/1000
ROH_file4$ROH_SUM <- as.numeric(ROH_file4$ROH_SUM)/1000000
ROH_file5$ROH_SUM <- as.numeric(ROH_file5$ROH_SUM)/1000000

# Add df together

Figure_S32 <- bind_rows(ROH_file2,ROH_file3,ROH_file4)
Figure_OSF <- bind_rows(ROH_file1,ROH_file5)

library(googlesheets4)  
# Google sheet is linked to Table S3 (Published ancients)
anc.info<-data.frame(read_sheet("https://docs.google.com/spreadsheets/d/1og33Y8H9PMYM2DTgDoM9vN63fRROa7v1E21twzl3UJc/edit?gid=1991801573#gid=1991801573",sheet = 3))

Figure_S32.merge<-left_join(Figure_S32, anc.info %>% select(Name, Country_continent_plot, Period_plot), by = c('Sample' = 'Name'),multiple = "any")
Figure_OSF.merge<-left_join(Figure_OSF, anc.info %>% select(Name, Country_continent_plot, Period_plot), by = c('Sample' = 'Name'),multiple = "any")

Figure_S32.merge$Country_continent_plot <- factor(Figure_S32.merge$Country_continent_plot,levels =c("Europe","Turkey","Israel","Armenia/Georgia","Iran","Turkmenistan/Uzbekistan"))
Figure_S32.merge$Sample <- factor(Figure_S32.merge$Sample,levels =c("Potterne1","Blagotin3","Blagotin2","Blagotin16","Blagotin1","Acem2","Acem1","Direkli6","Direkli1-2","Yoqneam2","Kazbegi1","Geor2","Hovk1","Azer4","Azer3-5","Ganjdareh18","Ganjdareh20","Ganjdareh22","Ganjdareh26","Ganjdareh3","Ganjdareh34","Abdul4","Lur12","Qazvin1","Fars4","Darre2","Semnan1-2","Semnan10","Semnan13","Semnan3","Semnan7","Semnan9","Monjukli8","Monjukli4","Bulak1","Bulak2"))
Figure_OSF.merge$Country_continent_plot <- factor(Figure_OSF.merge$Country_continent_plot,levels =c("Europe","Turkey","Israel","Armenia/Georgia","Iran","Turkmenistan/Uzbekistan"))
Figure_OSF.merge$Sample <- factor(Figure_OSD.merge$Sample,levels =c("Potterne1","Blagotin3","Blagotin2","Blagotin16","Blagotin1","Acem2","Acem1","Direkli6","Direkli1-2","Yoqneam2","Kazbegi1","Geor2","Hovk1","Azer4","Azer3-5","Ganjdareh18","Ganjdareh20","Ganjdareh22","Ganjdareh26","Ganjdareh3","Ganjdareh34","Abdul4","Lur12","Qazvin1","Fars4","Darre2","Semnan1-2","Semnan10","Semnan13","Semnan3","Semnan7","Semnan9","Monjukli8","Monjukli4","Bulak1","Bulak2"))


# Plot
colors <- c("#800000","#C21807","#FF817E","#D1EAF0","#0077B6","#012A4A")

plot_roh <- Figure_S32.merge %>%
  mutate(Bins=fct_relevel(Bins,'>16.0 Mb','8.0-16.0 Mb','4.0-8.0 Mb','2.0-4.0 Mb','1.0-2.0 Mb','0.5-1.0 Mb')) %>% 
  ggplot(aes(fill=Bins, x=Origin, y=ROH_SUM)) +
  geom_col(colour='black', size = 0.3,width = 0.55) +
  scale_fill_manual(values = colors, name="ROH bins") +
  labs(x= '', y = "ROH length (Mb)", size = 4) +
  facet_grid(Dataset ~ Country_continent_plot,scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0.1),limits = c(0,550)) +
  scale_x_discrete(expand = expand_scale(add = 0.3)) +
  theme_minimal() +
  theme(
    panel.margin = grid::unit(0.5, "lines"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major.y = element_line( size=.03, color="lightgray"),
    panel.grid.minor.y = element_line( size=.03, color="lightgray"),
    axis.text.y = element_text(size=4),
    axis.text.x = element_text(angle = 90,size=4),
    strip.text.x = element_text(size = 4),
    legend.key.width= unit(0.35, 'cm'),
    legend.key.height= unit(0.35, 'cm'),
    legend.text = element_text(size=5),
    legend.title = element_text(size=5),
    axis.title.y = element_text(size = 5)
  ) 

pdf("Figure_S32_ROH_comparison_individual_downsampled_plink_bcftools.pdf",width=8,height=4)
print(plot_roh)
dev.off()

# Plot OSF figure
plot_roh <- Figure_OSF.merge %>%
  mutate(Bins=fct_relevel(Bins,'>16.0 Mb','8.0-16.0 Mb','4.0-8.0 Mb','2.0-4.0 Mb','1.0-2.0 Mb','0.5-1.0 Mb')) %>% 
  ggplot(aes(fill=Bins, x=Origin, y=ROH_SUM)) +
  geom_col(colour='black', size = 0.3,width = 0.55) +
  scale_fill_manual(values = colors, name="ROH bins") +
  labs(x= '', y = "ROH length (Mb)", size = 4) +
  facet_grid(Dataset ~ Country_continent_plot,scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0.1),limits = c(0,550)) +
  scale_x_discrete(expand = expand_scale(add = 0.3)) +
  theme_minimal() +
  theme(
    panel.margin = grid::unit(0.5, "lines"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major.y = element_line( size=.03, color="lightgray"),
    panel.grid.minor.y = element_line( size=.03, color="lightgray"),
    axis.text.y = element_text(size=4),
    axis.text.x = element_text(angle = 90,size=4),
    strip.text.x = element_text(size = 4),
    legend.key.width= unit(0.35, 'cm'),
    legend.key.height= unit(0.35, 'cm'),
    legend.text = element_text(size=5),
    legend.title = element_text(size=5),
    axis.title.y = element_text(size = 5)
  ) 

pdf("Figure_S32_ROH_comparison_combined_dataset_plink_bcftools.pdf",width=8,height=4)
print(plot_roh)
dev.off()

# Load dataframes for Figure 5
ROH_file <- read.table("Figure5_ROH_sum_profiles_plink_bcftools_comparison.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

ROH_file$ROH_SUM <- as.numeric(ROH_file$ROH_SUM)/1000000

Figure_ROH.merge<-left_join(ROH_file, anc.info %>% select(Name, Country_continent_plot, Period_plot), by = c('Sample' = 'Name'),multiple = "any")

Figure_ROH.merge$Country_continent_plot <- factor(Figure_ROH.merge$Country_continent_plot,levels =c("Europe","Turkey","Israel","Armenia/Georgia","Iran","Turkmenistan/Uzbekistan"))
Figure_ROH.merge$Sample <- factor(Figure_ROH.merge$Sample,levels =c("Potterne1","Blagotin3","Blagotin2","Blagotin16","Blagotin1","Acem2","Acem1","Direkli6","Direkli1-2","Yoqneam2","Kazbegi1","Geor2","Hovk1","Azer4","Azer3-5","Ganjdareh18","Ganjdareh20","Ganjdareh22","Ganjdareh26","Ganjdareh3","Ganjdareh34","Abdul4","Lur12","Qazvin1","Fars4","Darre2","Semnan1-2","Semnan10","Semnan13","Semnan3","Semnan7","Semnan9","Monjukli8","Monjukli4","Bulak1","Bulak2"))

# Plot figure
plot_roh <- Figure_ROH.merge %>%
  mutate(Bins=fct_relevel(Bins,'>16.0 Mb','8.0-16.0 Mb','4.0-8.0 Mb','2.0-4.0 Mb','1.0-2.0 Mb','0.5-1.0 Mb')) %>% 
  ggplot(aes(fill=Bins, x=Origin, y=ROH_SUM)) +
  geom_col(colour='black', size = 0.3,width = 0.55) +
  scale_fill_manual(values = colors, name="ROH bins") +
  labs(x= '', y = "ROH length (Mb)", size = 4) +
  facet_grid(Dataset ~ Country_continent_plot,scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0.1),limits = c(0,550)) +
  scale_x_discrete(expand = expand_scale(add = 0.3)) +
  theme_minimal() +
  theme(
    panel.margin = grid::unit(0.5, "lines"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major.y = element_line( size=.03, color="lightgray"),
    panel.grid.minor.y = element_line( size=.03, color="lightgray"),
    axis.text.y = element_text(size=4),
    axis.text.x = element_text(angle = 90,size=4),
    strip.text.x = element_text(size = 4),
    legend.key.width= unit(0.35, 'cm'),
    legend.key.height= unit(0.35, 'cm'),
    legend.text = element_text(size=5),
    legend.title = element_text(size=5),
    axis.title.y = element_text(size = 5)
  ) 

pdf("Figure_S5_ROH_comparison_combined_and individual_dataset_plink_bcftools_combined.pdf",width=8,height=4)
print(plot_roh)
dev.off()
