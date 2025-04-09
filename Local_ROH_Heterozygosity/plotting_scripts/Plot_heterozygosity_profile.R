#!/usr/bin/env Rscript

# List required packages
list.of.packages <- c("ggplot2", "cowplot", "optparse")

# Install any missing packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# Load libraries
library(ggplot2)
library(cowplot)
library(optparse)

# Script to plot local heterozygosity profiles for given samples and chromosomes.
# Generates two profiles based on different plink/VCF heterozygosity sources.

# Set up command-line options
option_list = list(
  make_option(c("-s", "--sampleNames"), type = "character", default = NULL,
              help = "Sample name", metavar = "character"),
  make_option(c("-c", "--chr"), type = "character", default = NULL,
              help = "Chromosome number", metavar = "character")
)

# Parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Extract options
ind <- opt$sampleNames
chr <- opt$chr

### Load data files

# Load combined heterozygosity dataset
plink_HQ <- read.table(paste0("figure_hets_chr", chr, "_", ind, "_combined_glimpse2_GP99_MAF05_published_0-5X_geno0.vcf"),
                       header = FALSE, sep = ' ', comment.char = "")
plink_HQ1 <- read.table(paste0("hets_chr", chr, "_", ind, "_combined_glimpse2_GP99_MAF05_published_0-5X_geno0.vcf"),
                        header = FALSE, sep = ' ', comment.char = "")
						
# Load alternate 1.6M heterozygosity dataset
plink_1het <- read.table(paste0("figure_hets_chr", chr, "_", ind, "_glimpse2_GP99_MAF_0.05_TV_1.6M.vcf"),
                         header = FALSE, sep = ' ', comment.char = "")
plink_1het1 <- read.table(paste0("hets_chr", chr, "_", ind, "_glimpse2_GP99_MAF_0.05_TV_1.6M.vcf"),
                          header = FALSE, sep = ' ', comment.char = "")

### Plot 1: Combined HQ heterozygosity plot
plt1 <- ggplot() +
  geom_point(data = plink_HQ1, aes(x = V1, y = V3), color = 'gray', alpha = 0.9, shape = 21, fill = "gray", stroke = 0.2) +
  geom_point(data = plink_HQ, aes(x = V1, y = V3), color = 'black', shape = 21, fill = "black", stroke = 0.2) +
  labs(x = paste0('Chr', chr, ' (Bp)'), y = "Number of Heterozygotes per 50 SNPs", size = 8)

pdf(paste0("Het_profile_combined_chr", chr, "_", ind, ".pdf"), width = 4.5, height = 2.5)
plt1 + theme_classic() +
  theme(
    plot.margin = unit(rep(0, 4), "null"),
    axis.text.y = element_text(size = 8),
	panel.border = element_rect(colour='black', fill=NA, size=0.5),
    panel.grid = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0))
dev.off()

### Plot 2: 1.6M heterozygosity profile
plt2 <- ggplot() +
  geom_point(data = plink_1het1, aes(x = V1, y = V3), color = 'gray', alpha = 0.9, shape = 21, fill = "gray", stroke = 0.2) +
  geom_point(data = plink_1het, aes(x = V1, y = V3), color = 'black', shape = 21, fill = "black", stroke = 0.2) +
  labs(x = paste0('Chr', chr, ' (Bp)'), y = "Number of Heterozygotes per 50 SNPs", size = 8)

pdf(paste0("Het_profile_1.6M_chr", chr, "_", ind, ".pdf"), width = 4.5, height = 2.5)
plt2 + theme_classic() +
  theme(
    plot.margin = unit(rep(0, 4), "null"),
    axis.text.y = element_text(size = 8),
	panel.border = element_rect(colour='black', fill=NA, size=0.5),
    panel.grid = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0))
dev.off()
