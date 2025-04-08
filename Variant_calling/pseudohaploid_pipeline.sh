#!/bin/bash

# Pseudohaploid Pipeline  
# This script processes BAM files to generate pseudohaploid genotypes.  
# The pipeline includes:  
# 1. Listing BAM files that need to be pseudohaploidized  
# 2. Running a pseudohaploid calling script  
# 3. Performing quality checks on the generated haploid data  
# 4. Converting haploid data to PLINK format for downstream analysis  

# STEP 1: List BAM files that need to be pseudohaploidized  
ls /Ancients/BAM/*.bam | rev | cut -f 1 -d '/' | cut -f 2-50 -d '.' | rev > samples.txt  

# Iterate over each sample and process them  
while read i; do  
    # Create a filelist for each BAM file  
    ls /Ancients/BAM/"$i".bam > "$i".filelist  
    
    # STEP 2: Call pseudohaploid genotypes (can be parallelized)  
    /scripts/call_pseudohaploid_downsample.sh "$i"  
    
    # STEP 3: Perform quality control on haploid genotypes  
    python2 /scripts/SNP_ANGSD_LowCov_SanityCheck.py \
        /modern_reference_dataset/new_vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_angsd_pos.bim \
        "$i"_haploid.haplo.gz "$i"_haploid_filter.haplo.gz "$i"_problem_sites.txt $i  
    
    # Convert haploid data to PLINK format  
    /angsd/misc/haploToPlink "$i"_haploid_filter.haplo.gz "$i"_haploid_filter  
    
    # STEP 4: Adjust sample ID in PLINK .tfam file  
    # If the BAM filename matches the individual name, modify accordingly  
    sed -i s/ind0/$i_hap/g "$i"_haploid_filter.tfam  
    
    # STEP 5: Convert to PLINK binary format  
    /Software/Plink/plink --tfile "$i"_haploid_filter -chr-set 29 --make-bed --out "$i"_haploid_filter  
    
    # Clean up unnecessary files  
    rm "$i"_haploid_filter.tfam  

done < samples.txt  
