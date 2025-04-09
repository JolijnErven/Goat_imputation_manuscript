#!/bin/bash

# This script calculates the false positive rate (FPR) in 500kb genomic bins with
# 100kb steps, using two pre-binned BED files: one for heterozygous false positives
# and one for matched homozygous positions. It outputs a BED file containing:
#   - Chromosome
#   - Start
#   - End
#   - Total positions (false + true)
#   - False Positive Rate (%)

# Usage:
# ./calculate_false_positive_rate.sh <input_name> <chromosome> <sample_id> <coverage> <gp>

# Arguments:
#   $1 - Input file name prefix
#   $2 - Chromosome number
#   $3 - Sample ID
#   $4 - Coverage threshold
#   $5 - Genotype probability threshold


# Define input parameters
input_name=$1
chr=$2
sample=$3
cov=$4
gp=$5

#  Paste the two BED files side-by-side
# - Column 4 from file 1: False heterozygous counts
# - Column 8 from file 2: Homozygous match counts
#  If there are any false positives (>0), compute:
#   - Total calls (false + true)
#   - FPR = (false / total) * 100
paste \
    chr${chr}_${sample}_${cov}X_GP${gp}_${input_name}_False_heterozygous_positions_binned_500kb_100kb.bed \
    chr${chr}_${sample}_${cov}X_GP${gp}_${input_name}_Match_hom_binned_500kb_100kb.bed |
awk '{
    if ($4 > 0) {
        total = $4 + $8;
        fpr = ($4 / total) * 100;
        print $1, $2, $3, total, fpr;
    }
}' > chr${chr}_${sample}_${cov}X_GP${gp}_${input_name}_FPR_500kb_100kb.bed
