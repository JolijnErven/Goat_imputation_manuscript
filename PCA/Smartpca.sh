#!/bin/bash

# set parameters
dir="/modern_reference_dataset/"
dataset="new_vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans"

# Hardcoded for GP99, repeat for other GPs

# SmartPCA Pipeline for Genetic Analysis
# --------------------------------------
# This script performs the following steps:
# 1. Merges modern and ancient datasets.
# 2. Filters the merged dataset based on MAF and missingness criteria.
# 3. Converts PLINK format files to EIGENSTRAT format.
# 4. Adjusts EIGENSTRAT files.
# 5. Runs smartPCA for Principal Component Analysis (PCA).

# Set parameters
dir="/modern_reference_dataset/"
dataset="new_vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans"

# Merge modern dataset with ancient imputed dataset (ancients.txt - both downsampled imputed, downsampled pseudohaploid and high coverage) for GP99, repeat for other GPs
plink --noweb --chr-set 29 \
      --bfile ${dir}${dataset} \
      --merge-list ancients.txt \
      --recode \
      --out ${dataset}_ancient_merged_GP99

# Filter merged dataset for Minor Allele Frequency (MAF) >= 5% and missingness <= 1.5%
plink --file ${dataset}_ancient_merged \
      --chr-set 29 \
      --geno 0.015 \
      --recode \
      --extract ${dir}vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_MAF_0.05.pos.txt \
      --out ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015

# Convert PED format to EIGENSTRAT format
echo -e "genotypename: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.ped\n" \
        "snpname: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.map\n" \
        "indivname: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.ind\n" \
        "outputformat: EIGENSTRAT\n" \
        "genooutfilename: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.eigenstratgeno\n" \
        "snpoutfilename: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.snp\n" \
        "indoutfilename: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.ind\n" \
        "numchrom: 29\n" \
        "outliermode: 2" > ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.par

convertf -p ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.par \
         > conversion_${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.log

# Adjust EIGENSTRAT files
# Extract population IDs from individuals file
cut -f 1 -d'_' ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.ind | awk '{print $1}' > file.pop

# Modify individual file format for SmartPCA compatibility
awk '{print $1"\t"$2"\t"}' ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.ind \
    | paste - file.pop > ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015_alt.ind

# Remove populations/individuals that you want to project, which have to be removed from population file
grep -E -v "acem2|blagotin3|direkli1-2|semnan3" populations.txt | sort -u > populations_filtered.txt

# Create SmartPCA parameter file and run SmartPCA
echo -e "genotypename: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.eigenstratgeno\n" \
        "snpname: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.snp\n" \
        "indivname: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015_alt.ind\n" \
        "evecoutname: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015_projection.evec\n" \
        "evaloutname: ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.eval\n" \
        "poplistname: populations_filtered.txt\n" \
        "numoutlieriter: 0\n" \
        "killr2: YES\n" \
        "r2thresh: 0.2\n" \
        "lsqproject: YES\n" \
        "numchrom: 29\n" \
        "shrinkmode: YES\n" \
        "numthreads: 25" > ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015_projection.par

smartpca -p ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015_projection.par \
        >> projection_${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.015.log
