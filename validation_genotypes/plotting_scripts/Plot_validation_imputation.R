#!/usr/bin/env Rscript

# Load libraries
library(stringr)
library(dplyr)
library(MetBrewer)
require(ggplot2)


# Creates genotype validation plots (e.g. concordance and FPR) 
# Input files are in .csv format and are outputs from Concordance_FPR_FNR_calculations_MAF_threshold.r 
# and Concordance_FPR_FNR_calculations_MAF_trances.r scripts

# Set colours
pal<- met.brewer(name="Greek",n=5,type="discrete")

# ======================================================================================#
#                         MAF tranches-based plots										#
# ======================================================================================#

# Read dataframe
df <- data.frame(read.csv("correct_merged_validation_MAF_bins_concordance_type.csv"))

# Adjust strings
df <- df[df[,1] != "NAME", ]
df$TYPE <- as.character(df$TYPE)
df$TYPE <- gsub("trans", "Transition", df$TYPE)
df$TYPE <- gsub("tranvs", "Transversion", df$TYPE)
df$COVERAGE <- as.numeric(df$COVERAGE)
df$concord <- as.numeric(df$concord)
df$r2 <- as.numeric(df$r2)
df$MAF <- sapply(df$MAF, function(x) {
  values <- as.numeric(unlist(strsplit(x, "-")))
  paste0(values[1] * 100, "-", values[2] * 100)
})
df$MAF <- factor(df$MAF,levels = c("0-1","1-2","2-5","5-10","10-20","20-50"))

# Plot
pdf("MAF_bins_Het_r2_by_sample_GP_type_0-4x.pdf", width = 12,  height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=r2,group=interaction(DATASET,TYPE),colour=as.character(DATASET))) +
        geom_line(linewidth=1.0,aes(linetype=TYPE)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,4)) + facet_grid(NAME ~ MAF) +
        theme() + xlab("Coverage") + labs(subtitle = "MAF tranches %",color='Minimum GP',shape='Type') + ylab(expression("Non Reference R" ^ 2))  + scale_linetype_manual(values=c("solid","dotted")) +
theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

pdf("MAF_bins_Het_concordance_by_sample_GP_type_0-4x.pdf", width = 12,  height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=concord,group=interaction(DATASET,TYPE),colour=as.character(DATASET))) +
        geom_line(linewidth=1.0,aes(linetype=TYPE)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,4)) + facet_grid(NAME ~ MAF) +
        theme() + xlab("Coverage") + labs(subtitle = "MAF tranches %",color='Minimum GP',shape='Type') + ylab("Non Reference Concordance")  + scale_linetype_manual(values=c("solid","dotted")) +
theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

pdf("MAF_bins_Het_r2_by_sample_GP_type_0-1x.pdf", width = 12,  height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=r2,group=interaction(DATASET,TYPE),colour=as.character(DATASET))) +
  geom_line(linewidth=1.0,aes(linetype=TYPE)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,1)) + facet_grid(NAME ~ MAF) +
  theme() + xlab("Coverage") + labs(subtitle = "MAF tranches %",color='Minimum GP',shape='Type') + ylab(expression("Non Reference R" ^ 2))  + scale_linetype_manual(values=c("solid","dotted")) +
theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

pdf("MAF_bins_Het_concordance_by_sample_GP_type_0-1x.pdf", width = 12,  height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=concord,group=interaction(DATASET,TYPE),colour=as.character(DATASET))) +
  geom_line(linewidth=1.0,aes(linetype=TYPE)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,1)) + facet_grid(NAME ~ MAF) +
  theme() + xlab("Coverage") + labs(subtitle = "MAF tranches %",color='Minimum GP',shape='Type') + ylab("Non Reference Concordance")  + scale_linetype_manual(values=c("solid","dotted")) +
theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

# Repeat for Transitions and Transversions combined
df <- data.frame(read.csv("correct_merged_validation_MAF_bins_concordance.csv"))

df <- df[df[,1] != "NAME", ]
df$COVERAGE <- as.numeric(df$COVERAGE)
df$concord <- as.numeric(df$concord)
df$r2 <- as.numeric(df$r2)
df$MAF <- sapply(df$MAF, function(x) {
  values <- as.numeric(unlist(strsplit(x, "-")))
  paste0(values[1] * 100, "-", values[2] * 100)
})
df$MAF <- factor(df$MAF,levels = c("0-1","1-2","2-5","5-10","10-20","20-50"))

# Plot
pdf("MAF_bins_Heterozygous_+_Homo_alt_r2_by_sample_GP_0-4x.pdf", width = 12,  height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=r2,group=interaction(DATASET),colour=as.character(DATASET))) +
        geom_line(linewidth=1.0) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,4)) + facet_grid(NAME ~ MAF) +
        theme() + xlab("Coverage") + labs(subtitle = "MAF tranches %",color='Minimum GP',shape='Type') + ylab(expression("Non Reference R" ^ 2))  + scale_linetype_manual(values=c("solid","dotted"))
dev.off()

pdf("MAF_bins_Heterozygous_+_Homo_alt_r2_by_sample_GP_0-1x.pdf", width = 12,  height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=r2,group=interaction(DATASET),colour=as.character(DATASET))) +
        geom_line(linewidth=1.0) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,1)) + facet_grid(NAME ~ MAF) +
        theme() + xlab("Coverage") + labs(subtitle = "MAF tranches %",color='Minimum GP',shape='Type') + ylab(expression("Non Reference R" ^ 2))  + scale_linetype_manual(values=c("solid","dotted"))
dev.off()

pdf("MAF_bins_Heterozygous_+_Homo_alt_concordance_by_sample_GP_0-1x.pdf", width = 12,  height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=concord,group=interaction(DATASET),colour=as.character(DATASET))) +
        geom_line(linewidth=1.0) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,1)) + facet_grid(NAME ~ MAF) +
        theme(plot.subtitle = element_text(hjust = 0.5)) + xlab("Coverage") + labs(subtitle = "MAF tranches %",color='Minimum GP',shape='Type') + ylab("Non Reference Concordance")  + scale_linetype_manual(values=c("solid","dotted"))
dev.off()

pdf("MAF_bins_Heterozygous_+_Homo_alt_concordance_by_sample_GP_0-4x.pdf", width = 12,  height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=concord,group=interaction(DATASET),colour=as.character(DATASET))) +
        geom_line(linewidth=1.0) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,1)) + facet_grid(NAME ~ MAF) +
        theme(plot.subtitle = element_text(hjust = 0.5)) + xlab("Coverage") + labs(subtitle = "MAF tranches %",color='Minimum GP',shape='Type') + ylab("Non Reference Concordance")  + scale_linetype_manual(values=c("solid","dotted"))
dev.off()


# Read dataframes and plots for heterozygous FPR
df <- data.frame(read.csv("correct_merged_validation_MAF_bins_falseHetRate_type.csv"))

df <- df[df[,1] != "NAME", ]
df$TYPE <- as.character(df$TYPE)
df$TYPE <- gsub("trans", "Transition", df$TYPE)
df$TYPE <- gsub("tranvs", "Transversion", df$TYPE)
df$COVERAGE <- as.numeric(df$COVERAGE)
df$FPR <- as.numeric(df$FPR)
df$MAF <- sapply(df$MAF, function(x) {
  values <- as.numeric(unlist(strsplit(x, "-")))
  paste0(values[1] * 100, "-", values[2] * 100)
})
df$MAF <- factor(df$MAF,levels = c("0-1","1-2","2-5","5-10","10-20","20-50"))

pdf("MAF_bins_Heterozygous_FPR_GP_type_0-1X.pdf",width = 13, height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET))) +
        geom_point(data=subset(df, TYPE=="Transversion"),aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2)  +
        geom_point(data=subset(df, TYPE=="Transition"),aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2, stroke =0.4) + 
        scale_shape_manual(values = c(0,17)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0,1)) + facet_grid(NAME ~ MAF) + labs(subtitle = "MAF tranches %",x="Coverage",shape='Type', color='Minimum GP') + ylab("False-Positive Heterozygous rate (%)") +
        theme(plot.subtitle = element_text(hjust = 0.5)) 
dev.off()

pdf("MAF_bins_Heterozygous_FPR_GP_type_0-4X.pdf",width = 12, height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET))) +
  geom_point(data=subset(df, TYPE=="Transversion"),aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2) +
  geom_point(data=subset(df, TYPE=="Transition"),aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2, stroke =0.4) + 
  scale_shape_manual(values = c(0,17)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0,4)) + facet_grid(NAME ~ MAF) + labs(subtitle = "MAF tranches %", x="Coverage",color='Minimum GP',shape='Type') + ylab("False-Positive Heterozygous rate (%)") +
  theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

# Read dataframes and plots for heterozygous FNR
df <- data.frame(read.csv("correct_merged_validation_MAF_bins_false_negHetRate_type.csv"))

df <- df[df[,1] != "NAME", ]
df$TYPE <- as.character(df$TYPE)
df$TYPE <- gsub("trans", "Transition", df$TYPE)
df$TYPE <- gsub("tranvs", "Transversion", df$TYPE)
df$COVERAGE <- as.numeric(df$COVERAGE)

df$MAF <- sapply(df$MAF, function(x) {
  values <- as.numeric(unlist(strsplit(x, "-")))
  paste0(values[1] * 100, "-", values[2] * 100)
})
df$MAF <- factor(df$MAF,levels = c("0-1","1-2","2-5","5-10","10-20","20-50"))

pdf("MAF_bins_Heterozygous_FNR_GP_type_0-1X.pdf",width = 12, height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET))) +
  geom_point(data=subset(df, TYPE=="Transversion"),aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2)  +
  geom_point(data=subset(df, TYPE=="Transition"),aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2, stroke =0.4) + 
  scale_shape_manual(values = c(0,17)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0,1)) + facet_grid(NAME ~ MAF) + labs(subtitle = "MAF tranches %",x="Coverage",shape='Type', color='Minimum GP') + ylab("False-Negative Heterozygous rate (%)") +
  theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

pdf("MAF_bins_Heterozygous_FNR_GP_type_0-4X.pdf",width = 12, height = 7.32)
ggplot(data=df,aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET))) +
  geom_point(data=subset(df, TYPE=="Transversion"),aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2) +
  geom_point(data=subset(df, TYPE=="Transition"),aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2, stroke =0.4) + 
  scale_shape_manual(values = c(0,17)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0,4)) + facet_grid(NAME ~ MAF) + labs(subtitle = "MAF tranches %", x="Coverage",color='Minimum GP',shape='Type') + ylab("False-Negative Heterozygous rate (%)") +
  theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()


# ======================================================================================#
#                         MAF threshold-based plots										#
# ======================================================================================#

# Read dataframe
df <- data.frame(read.csv("header_final_merged_validation_MAF_threshold_concordance_type.csv"))


df <- df[df[,1] != "NAME", ]
df$TYPE <- as.character(df$TYPE)
df$TYPE <- gsub("trans", "Transition", df$TYPE)
df$TYPE <- gsub("tranvs", "Transversion", df$TYPE)
df$MAF <- as.character(df$MAF)
df$MAF <- gsub("0.05", "5% MAF threshold", df$MAF)
df$MAF <- gsub("0", "No MAF threshold", df$MAF)
df$MAF <- factor(df$MAF,levels = c("No MAF threshold","5% MAF threshold"))                                
df$COVERAGE <- as.numeric(df$COVERAGE)
df$concord <- as.numeric(df$concord)
df$r2 <- as.numeric(df$r2)


pdf("MAF_thresholds_het_r2_by_sample_GP_type_0-4x.pdf", width = 12,  height = 5.45)
ggplot(data=df,aes(x=COVERAGE,y=r2,group=interaction(DATASET,TYPE),colour=as.character(DATASET))) +
  geom_line(linewidth=1.0,aes(linetype=TYPE)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,4)) + facet_grid(MAF ~ NAME) +
  theme() + xlab("Coverage") + labs(color='Minimum GP',shape='Type') + ylab(expression("Non Reference R" ^ 2))  + scale_linetype_manual(values=c("solid","dotted")) 
dev.off()

pdf("MAF_thresholds_het_concordance_by_sample_GP_type_0-4x.pdf", width = 12,  height = 5.45)
ggplot(data=df,aes(x=COVERAGE,y=concord,group=interaction(DATASET,TYPE),colour=as.character(DATASET))) +
  geom_line(linewidth=1.0,aes(linetype=TYPE)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,4)) + facet_grid(MAF ~ NAME) +
  theme() + xlab("Coverage") + labs(color='Minimum GP',shape='Type') + ylab("Non Reference Concordance")  + scale_linetype_manual(values=c("solid","dotted")) 
dev.off()

pdf("MAF_thresholds_het_r2_by_sample_GP_type_0-1x.pdf", width = 12,  height = 5.45)
ggplot(data=df,aes(x=COVERAGE,y=r2,group=interaction(DATASET,TYPE),colour=as.character(DATASET))) +
  geom_line(linewidth=1.0,aes(linetype=TYPE)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,1)) + facet_grid(MAF ~ NAME) +
  theme() + xlab("Coverage") + labs(color='Minimum GP',shape='Type') + ylab(expression("Non Reference R" ^ 2))  + scale_linetype_manual(values=c("solid","dotted")) 
dev.off()

pdf("MAF_thresholds_het_concordance_by_sample_GP_type_0-1x.pdf", width = 12,  height = 5.45)
ggplot(data=df,aes(x=COVERAGE,y=concord,group=interaction(DATASET,TYPE),colour=as.character(DATASET))) +
  geom_line(linewidth=1.0,aes(linetype=TYPE)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0.1,1)) + facet_grid(MAF ~ NAME) +
  theme() + xlab("Coverage") + labs(color='Minimum GP',shape='Type') + ylab("Non Reference Concordance")  + scale_linetype_manual(values=c("solid","dotted")) 
dev.off()


# Read dataframes and plots for heterozygous FPR
df <- data.frame(read.csv("header_final_merged_validation_MAF_threshold_falseHetRate_type.csv"))

df$TYPE <- as.character(df$TYPE)
df$TYPE <- gsub("trans", "Transition", df$TYPE)
df$TYPE <- gsub("tranvs", "Transversion", df$TYPE)
df$MAF <- as.character(df$MAF)
df$MAF <- gsub("0.05", "5% MAF threshold", df$MAF)
df$MAF <- gsub("0", "No MAF threshold", df$MAF)
df$MAF <- factor(df$MAF,levels = c("No MAF threshold","5% MAF threshold"))
df$COVERAGE <- as.numeric(df$COVERAGE)
df$FPR <- as.numeric(df$FPR)

pdf("MAF_thresholds_Heterozygous_FPR_GP_type_0-1X.pdf",width = 12, height = 5.45)
ggplot(data=df,aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET))) +
  geom_point(data=subset(df, TYPE=="Transversion"),aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2)  +
  geom_point(data=subset(df, TYPE=="Transition"),aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2, stroke =0.6) + 
  scale_shape_manual(values = c(0,17)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0,1)) + facet_grid(MAF ~ NAME) + labs(x="Coverage",shape='Type', color='Minimum GP') + ylab("False-Positive Heterozygous rate (%)") 
dev.off()

pdf("MAF_thresholds_Heterozygous_FPR_GP_type_0-4X.pdf",width = 12, height = 5.45)
ggplot(data=df,aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET))) +
  geom_point(data=subset(df, TYPE=="Transversion"),aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2) +
  geom_point(data=subset(df, TYPE=="Transition"),aes(x=COVERAGE,y=FPR,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2, stroke =0.6) + 
  scale_shape_manual(values = c(0,17)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0,4)) + facet_grid(MAF ~ NAME) + labs(x="Coverage",color='Minimum GP',shape='Type') + ylab("False-Positive Heterozygous rate (%)") 
dev.off()

# Read dataframes and plots for heterozygous FNR
df <- data.frame(read.csv("header_final_merged_validation_MAF_threshold_false_negHetRate_type.csv"))

df$TYPE <- as.character(df$TYPE)
df$TYPE <- gsub("trans", "Transition", df$TYPE)
df$TYPE <- gsub("tranvs", "Transversion", df$TYPE)
df$MAF <- as.character(df$MAF)
df$MAF <- gsub("0.05", "5% MAF threshold", df$MAF)
df$MAF <- gsub("0", "No MAF threshold", df$MAF)
df$MAF <- factor(df$MAF,levels = c("No MAF threshold","5% MAF threshold"))
df$COVERAGE <- as.numeric(df$COVERAGE)

pdf("MAF_thresholds_Heterozygous_FNR_GP_type_0-1X.pdf",width = 12, height = 5.45)
ggplot(data=df,aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET))) +
  geom_point(data=subset(df, TYPE=="Transversion"),aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2)  +
  geom_point(data=subset(df, TYPE=="Transition"),aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2, stroke =0.6) + 
  scale_shape_manual(values = c(0,17)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0,1)) + facet_grid(MAF ~ NAME) + labs(x="Coverage",shape='Type', color='Minimum GP') + ylab("False-Negative Heterozygous rate (%)") 
dev.off()

pdf("MAF_thresholds_Heterozygous_FNR_GP_type_0-4X.pdf",width = 12, height = 5.45)
ggplot(data=df,aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET))) +
  geom_point(data=subset(df, TYPE=="Transversion"),aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2) +
  geom_point(data=subset(df, TYPE=="Transition"),aes(x=COVERAGE,y=falseHetRate,group=as.character(DATASET),colour=as.character(DATASET),shape=TYPE),size=2, stroke =0.6) + 
  scale_shape_manual(values = c(0,17)) + theme_bw(base_size=16)  + scale_colour_manual(values = pal) + xlim(c(0,4)) + facet_grid(MAF ~ NAME) + labs(x="Coverage",color='Minimum GP',shape='Type') + ylab("False-Negative Heterozygous rate (%)") 
dev.off()
