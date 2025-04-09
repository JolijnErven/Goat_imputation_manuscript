#!/usr/bin/env Rscript

# Load necessary packages
list.of.packages <- c("ggplot2", "cowplot", "optparse", "grid")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(ggplot2)
library(cowplot)
library(optparse)
library(grid)

# Script to plot local ROH profiles for given samples and chromosomes.

option_list = list(
    make_option(c("-s", "--sampleNames"), type="character", default=NULL,
                help="list of all sample names", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Base output", metavar="character"),
    make_option(c("-e", "--end"), type="character", default=NULL,
                help="Base output", metavar="character"),
    make_option(c("-c", "--chr"), type="character", default=NULL,
                help="Base output", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

ind <- opt$sampleNames
chr <- opt$chr
chr_end <- opt$end

### File reading
# plink
plink <- read.table("final_size_combined_glimpse2_GP99_MAF05_published_0-5X_geno0_500_100kb_1het_200SNP.hom", header=T, sep='\t',comment.char = "")
plink_1.2M <- read.table(paste0("final_size_", ind, "_glimpse2_GP99_MAF_0.05_TV_1.2M_500kb_100kb_200SNPs_1het.hom"), header=T, sep='\t',comment.char = "")
plink_1.6M <- read.table(paste0("final_size_", ind, "_glimpse2_GP99_MAF_0.05_TV_1.6M_500kb_100kb_200SNPs_1het.hom"), header=T, sep='\t',comment.char = "")

# Get correct chromosome
plink.chr <- plink[plink$CHR == chr,]
plink_1.2M.chr <- plink_1.2M[plink_1.2M$CHR == chr,]
plink_1.6M.chr <- plink_1.6M[plink_1.6M$CHR == chr,]

# Get sample
plink.ind <- plink.chr[plink.chr$IID == ind,]

# Bcftools
b1 <- read.table("final_size_combined_published_0-5X_10Q_200SNPs.txt", header=T, sep='\t',comment.char = "")
b_1.6M <- read.table(paste0("final_size_", ind, "_glimpse2_GP99_MAF_0.05_TV_1.6M_Q10_200SNPs.txt"), header=T, sep='\t',comment.char = "")

# Get correct chromosome
b1.chr <- b1[b1$CHR == chr,]
b_1.6M.chr <- b_1.6M[b_1.6M$CHR == chr,]

# Get sample
b1.ind <- b1.chr[b1.chr$ID == ind,]

# Get ROHan
rohan  <- read.table(paste0("ROHan_", ind, ".txt"), header=T, sep='\t',comment.char = "")
rohan.chr <- rohan[rohan$CHROM == chr,]


# Plot ROH segments -- add more files if necessary 
plt <- ggplot() +
  geom_segment(data=rohan.chr, mapping=aes(x=BEGIN, y=0.0000032, xend=END, yend=0.0000032), color=rohan.chr$X, size=2, lineend="butt") +
  geom_segment(data=plink.ind, mapping=aes(x=as.numeric(POS1), y=0.0000025, xend=as.numeric(POS2), yend=0.0000025), color=plink.ind$X, size=2, lineend="butt") +
  geom_segment(data=plink_1.2M.chr, mapping=aes(x=POS1, y=0.000002, xend=POS2, yend=0.000002), color=plink_1.2M.chr$X, size=2, lineend="butt") +
  geom_segment(data=plink_1.2M.chr, mapping=aes(x=as.numeric(POS1), y=0.0000015, xend=as.numeric(POS2), yend=0.0000015), color=plink_1.2M.chr$X, size=2, lineend="butt") +
  geom_segment(data=b1.ind, mapping=aes(x=as.numeric(Start), y=0.00000075, xend=as.numeric(End), yend=0.00000075), color=b1.ind$X, size=2, lineend="butt") +
  geom_segment(data=b_1.6M.chr, mapping=aes(x=Start, y=0.00000025, xend=End, yend=0.00000025), color=b_1.6M.chr$X, size=2, lineend="butt") +
  labs(x= paste0('Chr', chr, ' (Bp)'), y = "", size = 8) +
  geom_segment(data=plink.ind, mapping=aes(x=0, y=0.00000000002, xend=as.numeric(chr_end), yend=0.000000000002), color="white", size=0.00038, lineend="butt") 

# Save plot to PDF
pdf(paste0("ROH_profiles_chr", chr, "_", ind, ".pdf"), width = 5, height = 2.5)
plt + 
  theme_classic() + 
  theme(
    plot.margin = unit(rep(0, 4), "null"),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5) 
  ) + 
  scale_y_continuous(
    limits = c(0, 0.0000035), expand = c(0, 0),
    breaks = c(0.00000025, 0.00000075, 0.0000015, 0.000002, 0.0000025, 0.0000032),
    labels = c("Bcftools individual 1.6M TV", "Bcftools >0.5x combined", 
               "Plink individual 1.6M TV", "Plink individual 1.2M TV", 
               "Plink >0.5x combined", "ROHan")
  ) + 
  scale_x_continuous(limits = c(0, as.numeric(chr_end)), expand = c(0, 0))

dev.off()
