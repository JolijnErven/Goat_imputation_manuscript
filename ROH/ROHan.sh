#!/bin/bash

# ROHan Pipeline for Detecting Runs of Homozygosity
# --------------------------------------------------
# This script performs the following steps:
# 1. Runs samtools to apply base quality score recalibration.
# 2. Indexes the BAM file for faster processing.
# 3. Estimates DNA damage using ROHan's estimateDamage.pl.
# 4. Runs ROHan to detect runs of homozygosity.

# Input: BAM file prefix (e.g., sample name)
dir="/goat/ancients_ARS1/"
input_prefix=$1

# Step 1: Recalibrate base quality scores
samtools calmd -@ 20 -b ${dir}${input_prefix}.bam \
    /Reference_Genomes/goat/ARS1.fa > ${input_prefix}.bam

# Step 2: Index the BAM file
samtools index -@ 20 ${input_prefix}.bam

# Step 3: Estimate DNA damage patterns
/Software/rohan/src/estimateDamage.pl \
    /Reference_Genomes/goat/ARS1.fa ${input_prefix}.bam

# Step 4: Run ROHan to detect runs of homozygosity
/raid_md0/Software/rohan/src/rohan \
    --rohmu 5.7e-5 -t 30 --auto chrom.txt \
    --deam5p ${input_prefix}.5p.prof \
    --deam3p ${input_prefix}.3p.prof \
    -o ${input_prefix} \
    /Reference_Genomes/goat/ARS1.fa ${input_prefix}.bam
