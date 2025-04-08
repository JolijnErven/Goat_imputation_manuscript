#!/bin/bash

# Script to downsample BAM files to specific coverages.
# This script iterates over a predefined set of downsampling fractions and 
# uses Samtools to generate BAM files with reduced read coverage.
# Usage: ./script.sh <input_bam_basename> (without .bam extension)

# Check if $1_stats directory exists. If it does, skip running bamqc.
if [ ! -d "$1_stats/genome_results.txt" ]; then
    echo "Running Qualimap BAMQC for $1.bam"
    /Software/qualimap_v2.2.1/qualimap bamqc --bam $1.bam --java-mem-size=5G -nt 9 --outdir $1_stats
else
    echo "Stats already exist for $1, skipping Qualimap BAMQC."
fi

# Extract mean coverage from the genome_results.txt file
COV=$(grep "mean coverageData" $1_stats/genome_results.txt | cut -f 2 -d "=" | sed 's/X//g')

# Loop through predefined downsampling fractions, adjust if necesarry
for i in {0.1,0.15,0.2,0.25,0.3,0.40,0.50,0.75,1,2,4}; do 
    # Calculate the downsampling fraction
    frac=$((i/COV))
    
    echo "Downsampling to $i with fraction $frac for $1"
    # Use Samtools to downsample the BAM file
    /Software/samtools-1.13/samtools view -s $frac $1.bam -o $1_"$i"X.bam 
done
