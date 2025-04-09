#!/bin/bash

# This script extracts SNP positions within a specified genomic region
# from a reference SNP position file

# Arguments:
#   $1 - Chromosome (e.g., "2")
#   $2 - Start position of region (e.g., "1000000")
#   $3 - End position of region (e.g., "2000000")


# Define input parameters
chr=$1
start=$2
end=$3

# Filter SNPs by chromosome and position range, and output in chr,pos format
awk -v chr=$chr -v start=$start -v end=$end -v OFS='\t' '
    $1 == chr && $2 >= start && $2 <= end {
        print $1 "," $2
    }
' /modern_reference_dataset/vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans.pos.txt
