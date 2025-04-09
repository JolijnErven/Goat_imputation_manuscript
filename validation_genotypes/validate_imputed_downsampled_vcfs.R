#!/usr/bin/env Rscript

# Made to prep imputed VCF files for concordance calculations
# Created by Kevin Daly -- updated by Jolijn Erven

# Load required libraries
require(stringr)
library(dplyr)

# Read command line arguments
args = commandArgs(trailingOnly=TRUE)

# Calling of script is like this (adjust for your own files);
# for CHR in {1..29}; do for TRUE in $(ls path_to_TRUTH_VCF/*chr${CHR}_*.vcf.gz); 
# do SAMPLE=`echo $TRUE | rev | cut -f1 -d'/' | rev | cut -f1 -d'_'`; for DOWNSAMPLED in $(ls ${SAMPLE}*chr${CHR}_*vcf| cut -f1 -d'.'); do SAMPLE=`echo $DOWNSAMPLED | 
# cut -f1 -d'_'`; COV=`echo $DOWNSAMPLED | cut -f2 -d'_' | cut -f1 -d'.' | sed -e "s/-/./g" | sed -e "s/X//g"` ; GP=`echo $DOWNSAMPLED | cut -f3 -d'_' | cut -f1 -d'.'` ;echo Rscript THIS_SCRIPT.r $TRUE ${DOWNSAMPLED}.vcf ${SAMPLE}
# $COV ${DOWNSAMPLED}.validation "2>" ${DOWNSAMPLED}.error >> parallel_validation.sh

# args1 is the TRUTH high coverage VCF
# args2 is the imputed downsampled VCF
# args3 is sample name
# args4 is coverage
# args5 is genotype probability (GP) threshold in imputed VCF
# args6 is outputfile name

# Read input files
vcf.true<-read.table(args[1])
vcf.sample<-read.table(args[2])

# Extract genotype information (first field from column 10) for the datasets
vcf.true$V11<-str_split_fixed(vcf.true$V10, pattern = ":", n=2)[,1]
vcf.sample$V11<-str_split_fixed(vcf.sample$V10, pattern = ":", n=2)[,1]


# Extract genotype probability (GP) from the INFO field (column 8) of the sample dataset
vcf.sample$V12<-as.numeric(str_split_fixed(str_split_fixed(vcf.sample$V8, pattern = ";", n=2)[,1],"=",n=2)[,2])

# Convert data frames
vcf.sample<-data.frame(vcf.sample)
vcf.true<-data.frame(vcf.true)

# Merge true and sample data based on chromosome and position (columns V1 and V2)
vcf.merged<-merge(vcf.sample, vcf.true, by=c("V1","V2"))[c("V1","V2","V4.x","V5.x","V12","V11.x","V11.y")]

# Add sample name from command-line arguments
vcf.merged$sample<-args[3]
vcf.merged$coverage<-args[4]

# Process coverage information
vcf.merged$coverage<- str_replace(vcf.merged$coverage, "-",".")
vcf.merged$coverage<-as.numeric(vcf.merged$coverage)

# Add genotype probability (GP) threshold from command-line arguments
vcf.merged$GP<-as.numeric(args[5])

# Write the processed data to an output file
write.table(x = vcf.merged,file=paste0(args[6]),row.names = F, col.names = F, sep = ",")

