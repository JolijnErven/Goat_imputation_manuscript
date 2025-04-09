#!/usr/bin/env Rscript

# Load libraries
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)

# This script generates Manhattan plots from genomic BED files representing 
# local concordance measures. Input files are output files from Local_concordance_calculations.sh
# and Local_heterozygous_FPR_calculations.sh, files are individual bed files with prefix Merged_
# So not the combined
# Two types of files are processed:
# 1. Non-reference local concordance data 
# 2. Heterozygous false-positive rate (FPR) concordance data.

mypalette <- c("#00ced1","black")

gg.manhattan <- function(df, threshold, col, ylims, title,ytitle){
  df.tmp <- df %>% 
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    left_join(df, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot)
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=FST)) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=1, size=0.25,stroke=0.05) +
    scale_color_manual(values = rep(col, 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome", y=paste0(ytitle)) +
    
    # add genome-wide sig lines
    geom_hline(yintercept = sig,linetype="dashed",size=0.2) +

    theme_bw(base_size = 4) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.grid.major.y = element_line(colour='gray', size =0.1),
      panel.grid.minor.y = element_line(colour='gray', size =0.1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

# Loop through Non Reference local concordance bed files
files = read.table('get_local_het_homoalt_conc.txt')
listfiles = as.list(files)

for (x in files){
  print(x)
  for (id in x){
    df = read.table(paste0(id, '.bed'), header = TRUE, sep=" ")
    colnames(df)<-c("CHR","BP","END","COUNT","FST")
    sig =  quantile(df$FST, c(0.01))
    
    tiff(paste0(id, '.tiff'),units="in", width=4, height=1.3, res=1900)
    print(gg.manhattan(df, threshold=sig, col=mypalette, ylims=c(0,1), title="", ytitle='Non Reference Concordance'))
    dev.off()
    
    pdf(paste0(id, '.pdf'), width=4, height=1.3)
    print(gg.manhattan(df, threshold=sig, col=mypalette, ylims=c(0,1),title=""))
    dev.off()
  }
}

# Loop through heterozygous FPR local concordance bed files
files = read.table('get_local_FPR.txt')
listfiles = as.list(files)

for (x in files){
  print(x)
  for (id in x){
    df = read.table(paste0(id, '.bed'), header = TRUE, sep=" ")
    colnames(df)<-c("CHR","BP","END","COUNT","FST")
    sig =  quantile(df$FST, c(0.99))
    
    tiff(paste0(id, '.tiff'),units="in", width=4, height=1.3, res=1900)
    print(gg.manhattan(df, threshold=sig, col=mypalette, ylims=c(0,1), title="", ytitle='Heterozygous False-Positive rate (%)'))
    dev.off()
    
    pdf(paste0(id, '.pdf'), width=4, height=1.3)
    print(gg.manhattan(df, threshold=sig, col=mypalette, ylims=c(0,100),title=""))
    dev.off()
  }
}
