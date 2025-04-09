#!/bin/bash

# Script to perform local concordance analysis in sliding windows for chromosomes 1-29
# This script consists of three main steps:
# 1. Divide sites into sliding windows
# 2. Calculate concordance for each window
# 3. Merge results and extract low concordance regions

# STEP 1: Divide sites into sliding windows
# - run_bed.txt consists of: Sample ID, Coverage, GP
# - make_bed.txt consists of: Prefix of .csv position file (generated from Concordance_FPR_FNR_calculations_MAF_threshold_local_concordance.r)

while read i; do 
    while read d; do 
        for c in {1..29}; do 
            echo /scripts/get_bed_windows.sh "$i" "$c" "$d" 
        done 
    done < run_bed.txt 
done < make_bed.txt > new_run_parallel_bed.sh

# Run bed window generation in parallel (20 concurrent jobs)
parallel -a new_run_parallel_bed.sh -j20

# STEP 2: Calculate concordance for every window
# - run_bed.txt consists of: Sample ID, Coverage, GP
# - make_bed.txt consists of: Prefix of .csv position file (from Concordance_FPR_FNR_calculations_MAF_threshold_local_concordance.r)

while read i; do 
    while read d; do 
        for c in {1..29}; do 
            echo "$i" "$c" "$d" 
        done 
    done < run_bed.txt 
done < make_bed.txt \
| sed 's/_Match_/ /g' \
| sed 's/_Discordant_/ /g' \
| sort | uniq > new_get_calc_bins.txt

# Generate parallel execution commands for concordance calculations
while read i; do 
    echo /scripts/calc_FPR_bins.sh "$i"
done < new_get_calc_bins.txt | grep -v positions > new_run_parallel_bins.sh

# Run concordance calculations in parallel (10 concurrent jobs)
parallel -a new_run_parallel_bins.sh -j10

# STEP 3: Merge results from all chromosomes for each sample
# Extract relevant parts of filenames (from chr9 as a reference)
ls *chr9*_FPR_500kb_100kb.bed | cut -f 2-50 -d "_" > concat.txt

# Concatenate results from all chromosomes into one merged file per sample
while read i; do 
    cat chr{1..29}_"$i" > Merged_"$i"
done < concat.txt

# STEP 4: Filter merged bed file for 1 percentile and merge adjoining regions
while read i; do
    # Get number of lines when looking at 1 percentile
    perc=$(wc -l Merged_"$i".bed | awk '{print int($1 * 0.01)}')
    # Get lowest 1 percentile of concordance regions ; note -r in sort because we want the highest value not the lowest
    sort -k5 -n -r Merged_"$i".bed | head -${perc} | sort -k 1,1 -k 2,2 -n | sed 's/ /\t/g' > Merged_"$i"_99.9_perc.bed
    # Merge adjacent regions
    bedtools merge -i Merged_"$i"_99.9_perc.bed > Merged_"$i"_99.9_perc_region.bed

done < concat.txt

# STEP 5: Get low concordance regions in the majority of the samples and add gene information
# Get files to iterate over removing sample specific information; change sample IDs accordingly
ls Merged_*FPR_500kb_100kb_99.9_perc_region.bed | cut -f 3-50 -d "_" > get_majority.txt

# Combine test individuals and retrieve low concordance regions in the majority of individuals
while read i; do
    # Merge test individuals and keep count of presence of low concordance regions
    bedtools multiinter -names Acem Blagotin Direkli Semnan -header -i Merged_acem2_"$i" Merged_blagotin3_"$i" Merged_direkli1-2_"$i" Merged_semnan3_"$i" > combined_"$i"
    # Filter for regions where the majority (>2) have low concordance regions
    awk '{if ($4>2) print $0}' combined_"$i" > majority_"$i"
    
    # Get genes within and closest to the windows 
    bedtools closest -t all -a majority_"$i" -b new_GCF_001704415.1_ARS1_genes.bed -g /Reference_Genomes/goat/ARS1.fa.fai | cut -f 1-5,13 | uniq \
    | bedtools merge -c 4,5,6 -o distinct,distinct,distinct > majority_genes_"$i"

done < get_majority.txt

# Afterwards, genes were visually inspected
