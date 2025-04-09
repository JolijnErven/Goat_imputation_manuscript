#!/bin/bash

# Description:
# This script processes a CSV file containing genomic data for a specific sample
# and generates a binned coverage file using bedtools. 
# The output is a `.bed` file with coverage mapped to 500kb windows,
# sliding every 100kb.

# Usage:
# ./bin_genomic_coverage.sh <input_name> <chromosome> <sample_id> <coverage> <gp>

# Arguments:
#   $1 - Input file name prefix (e.g., "run1" if input file is "chr2.run1.csv")
#   $2 - Chromosome number (e.g., "2")
#   $3 - Sample ID
#   $4 - Coverage threshold
#   $5 - Genotype probability threshold


# Define input parameters
input_name=$1
chr=$2
sample=$3
cov=$4
gp=$5

# Transform CSV to bed.
awk -F',' -v OFS='\t' -v sample=$sample -v cov=$cov -v gp=$gp -v chr=$chr '
    $8 == sample && $1 == chr && $9 == cov && $10 == gp {
        print $1, $2, $2+1, 1
    }
' chr${chr}.${input_name}.csv | sort -k2,2n > chr${chr}_${sample}_${cov}X_GP${gp}_${input_name}.bed

# Generate genomic windows (500kb size, 100kb step) and map coverage counts
# to these windows. Replace any missing values (".") with 0.
bedtools makewindows \
    -g /raid_md0/jolijn/reference/chr${chr}_Goat_chrom.txt \
    -w 500000 -s 100000 |
bedtools map \
    -a - \
    -b chr${chr}_${sample}_${cov}X_GP${gp}_${input_name}.bed \
    -c 4 -o sum |
sed 's/\./0/g' > chr${chr}_${sample}_${cov}X_GP${gp}_${input_name}_binned_500kb_100kb.bed
