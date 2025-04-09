#!/bin/bash

# ROH Exploration with PLINK for 0.5X Imputed all tested individuals together
# The analysis produced Figure S30 in the study.

# Merge and filter high coverage sample with imputed sample (GP99 - 0.5X)
ls /Imputed/glimpse_Plink/*0.5X*GP99.bed | sed 's/.bed//g' > merge_list.txt
ls /high_coverage/*_Q30_HC.bed | sed 's/.bed//g' | grep -v acem2_Q30_HC >> merge_list.txt

mkdir bcftools_ROH_no_missingness 

plink --bfile /high_coverage/acem2_Q30_HC --chr-set 29 --geno 0 --merge-list merge_list.txt --recode vcf
--extract /modern_reference_dataset/vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_MAF_0.05.pos.txt
--out bcftools_ROH_no_missingness/Combined_0.5X_HC_MAF_0.05_no_missingness

cd bcftools_ROH_no_missingness 
bgzip Combined_0.5X_HC_MAF_0.05_no_missingness.vcf
tabix Combined_0.5X_HC_MAF_0.05_no_missingness.vcf.gz

# Fix alleles in vcf (Due to plink switching alleles)
/Software/bcftools-1.17/bcftools norm -f Reference_Genomes/goat/ARS1.fa \
    --check-ref ws --do-not-normalize -Oz -o Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.vcf.gz \
	Combined_0.5X_HC_MAF_0.05_no_missingness.vcf.gz

zcat Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.vcf.gz | grep "#" > Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix.vcf
zcat Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.vcf.gz | awk '$4 != $5' | grep -v "#" >> Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix.vcf
bgzip Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix.vcf
tabix Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix.vcf.gz

# Calculate ROHs with bcftools
/Software/bcftools-1.17/bcftools roh -G 30 --AF-dflt 0.4 Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix.vcf.gz \
-o Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix.txt

# Filter ROH for minmum 200SNPs, 10 quality and a length of 500kb 
grep RG Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix.txt | awk '{if ($8>10) print $0}' |  awk '{if ($7>=200) print $0}' \
 | awk '{if ($6>=500000) print $0}' > Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q.txt

# Transform .hom file into ROH bin sizes and ease for plotting
python3 /scripts/get_ROH_demographics_BCFtools.py -i Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q.txt \
-o Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q
python3 /scripts/Transform_ROH_demographics.py -i ROH_demo_plots_Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q.txt \
-o ROH_demo_plots_Combined_0.5X_HC_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q