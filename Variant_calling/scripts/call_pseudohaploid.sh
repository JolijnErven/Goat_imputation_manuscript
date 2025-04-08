#!/bin/bash

# Script to perform haplotype calling using ANGSD  
# This script takes a BAM file list as input and generates haploid genotype calls.  
# The output includes allele counts, filtered by quality and mapping criteria.  

# Check if binary sites .txt.bin file is present  
if [ ! -f /modern_reference_dataset/new_vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_angsd_pos.txt.bin ]; then
    echo "Binary sites file not found. Generating it now..."
    /Software/angsd/angsd sites index /modern_reference_dataset/new_vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_angsd_pos.txt 
else
    echo "Binary sites file found. Proceeding with analysis..."
fi

# Command to run ANGSD for haplotype calling  
/Software/angsd/angsd \  
    -doHaploCall 1 \          # Perform haplotype calling  
    -doCounts 1 \             # Count allele occurrences  
    -dumpCounts 1 \           # Output allele count data  
    -C 50 \                   # Adjust BAQ (Base Alignment Quality) threshold  
    -bam $1.filelist \        # Input BAM file list  
    -out $1_haploid \         # Output prefix  
    -uniqueOnly 1 \           # Consider only uniquely mapped reads  
    -trim 4 \                 # Trim bases from both ends of reads  
    -doMajorMinor 3 \         # Infer major/minor alleles  
    -remove_bads 1 \          # Remove bad reads  
    -minMapQ 30 \             # Minimum mapping quality threshold  
    -MinQ 25 \                # Minimum base quality threshold  
    -rf chrom.txt \           # Restrict to specific chromosomes  
    -nThreads 10 \            # Number of threads to use  
    -sites /modern_reference_dataset/new_vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_angsd_pos.txt \  # List of variant sites of the reference panel used for imputation
    -ref /Reference_Genomes/goat/ARS1.fa  # Reference genome  
