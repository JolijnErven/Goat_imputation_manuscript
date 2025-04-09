#!/bin/bash

# ROH Exploration with PLINK for 0.25X,0.5X,0.75X,1X Imputed individual specific datasets
# The analysis produced Figure S29 in the study.

# Merge and filter high coverage sample with imputed samples (GP99 - 0.25X,0.5X,0.75X,1X)
echo -e 'acem1\nblagotin3\ndirekli1-2\nsemnan3' > samples.txt
# Merge and filter for MAF 0.05 and no missingness

mkdir -p ROH_no_missingness

while read i; do
   ls /Imputed/glimpse_Plink/"$i"_GP99.bed | grep -E "0.25X|0.5X|0.75X|1X" | sed 's/.bed//g' > merge_list.txt
   plink --bfile /high_coverage/"$i"_Q30_HC --chr-set 29 --geno 0 --merge-list merge_list.txt --make-bed
   --extract /modern_reference_dataset/vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_MAF_0.05.pos.txt
   --out ROH_no_missingness/"$i"_MAF_0.05_no_missingness
; done 

# Calculate ROH with a minimum of 200 SNPs
while read i; do
    plink --bfile ROH_no_missingness/"$i"_MAF_0.05_no_missingness \
    --chr-set 29 \
    --homozyg \
    --homozyg-density 50 \
    --homozyg-gap 100 \
    --homozyg-kb 500 \
    --homozyg-snp 200 \
    --homozyg-window-het 1 \
    --homozyg-window-snp 50 \
    --homozyg-window-threshold 0.05 \
    --memory 40000 \
    --threads 4 \
    --out ROH_no_missingness/"$i"_MAF_0.05_no_missingness_500kb_100kb_1het_200SNPs
; done < samples.txt

# Get ROH Bin statistics and transpose
cd ROH_no_missingness
ls *.hom | sed 's/.hom//g' > hom.samples.txt

while read i; do 
    python3 /scripts/get_ROH_demographics.py -i "$i".hom -o $i
	
	python3 /scripts/Transform_ROH_demographics.py -i ROH_demo_plots_"$i".txt -o ROH_demo_plots_"$i"

; done < hom.samples.txt

# Prepare for plotting
# Cat all samples
cat *sum_transpose.txt | grep -v Sample > Combined_samples_no_missingness_ROH_500kb_100kb_1het_200SNPs.txt
# Make sample and coverage a seperate column
sed -i 's/_/\t/g' Combined_samples_no_missingness_ROH_500kb_100kb_1het_200SNPs.txt

# Add header
sed -i '1s/^/Sample\tCoverage\tBins\tROH_SUM\n/' Combined_samples_no_missingness_ROH_500kb_100kb_1het_200SNPs.txt
