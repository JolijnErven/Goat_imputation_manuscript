#!/bin/bash


# Description:
# This script calculates the non-reference concordance rate in 500kb genomic bins
# (sliding by 100kb) using two pre-binned BED files, One for matched non-reference genotypes
#   and one for discordant non-reference genotypes. It outputs a BED file containing:
#   - Chromosome
#   - Start
#   - End
#   - Total positions (match + discordant)
#   - Concordance rate (match / total)

# Usage:
# ./calculate_nonref_concordance.sh <input_name> <match_type> <chromosome> <sample_id> <coverage> <gp>

# Arguments:
#   $1 - Input file name prefix
#   $2 - Match type (e.g., "Het_hom_alt" or other identifier)
#   $3 - Chromosome number
#   $4 - Sample ID
#   $5 - Coverage threshold
#   $6 - Genotype probability threshold


# Define input parameters
input_name=$1
match_type=$2
chr=$3
sample=$4
cov=$5
gp=$6

# Step 1: Paste the matched and discordant BED files side-by-side
# Step 2: If there are any matched non-reference positions (>0), calculate:
#   - Total calls (matched + discordant)
#   - Concordance rate = matched / total
paste \
    chr${chr}_${sample}_${cov}X_GP${gp}_${input_name}_Match_${match_type}_binned_500kb_100kb.bed \
    chr${chr}_${sample}_${cov}X_GP${gp}_${input_name}_Discordant_${match_type}_binned_500kb_100kb.bed |
awk '{
    if ($4 > 0) {
        total = $4 + $8;
        concordance = $4 / total;
        print $1, $2, $3, total, concordance;
    }
}' > chr${chr}_${sample}_${cov}X_GP${gp}_${input_name}_${match_type}_concordance_500kb_100kb.bed
