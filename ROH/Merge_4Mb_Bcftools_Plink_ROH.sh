#!/bin/bash

# Script to merge ROH profiles from Plink and Bcftools
# Script merges large (>4Mb) ROH from Bcftools with ROHs from Plink


# Create directory
mkdir -p plink_bcftools

# Get all individual 1.6M bcftools ROH files
ls *_MAF_0.05_TV_1.6M_200SNPs_10Q.txt | grep -v demo  | sed 's/.txt//g' > files.txt

# Filter bcftools ROH for only ROH >=4Mb and transform into bed format
while read i; do awk -v OFS='\t' '{ if ($6>=4000000) print $3,$4,$5,$2,$6}' "$i".txt > ../plink_bcftools/"$i"_4mb.bed ; done < files.txt

cd plink_bcftools

# Transform plink 1.6M .hom files to .bed files
ls /published_paleogenomes_ROH/*_MAF_0.05_TV_1.6M_500kb_100kb_200SNPs_1het.hom | rev| cut -f 1 -d '/' | rev  | sed 's/.hom//g' > files.txt

while read i; do awk -v OFS='\t' '{print $4,$7,$8,$1,$9}' /published_paleogenomes_ROH/"$i".hom > "$i".bed ; done < files.txt

# Merge ROH segments and keep information on how many and size of ROHs merged
# Get sample IDs, Sample ID is in the first segment of file name
ls *MAF_0.05_TV_1.6M_500kb_100kb_200SNPs_1het.bed | cut -f 1 -d "_" > samples.txt

while read i; do cat "$i"_glimpse2_GP99_MAF_0.05_TV_1.6M_Q10_200SNPs_4mb.bed "$i"_glimpse2_GP99_MAF_0.05_TV_1M_500kb_100kb_200SNPs_1het.bed | grep -v CHR | \
awk -v OFS='\t' '{print $1,$2,$3,$4,$5}' | sort -k1,1 -k2,2 -n | bedtools merge -c 4,5 -o count,collapse | awk -v OFS='\t' -v id=${i} '{print id,$1,$2,$3,$4,$5,$3-$2}' \
 > "$i"_GP99_MAF_0.05_TV_1.6M_bcftools_plink_ROH.bed ; done < samples.txt
 
# calculate bin size + transform
ls *_GP99_MAF_0.05_TV_1.6M_bcftools_plink_ROH.bed |sed 's/.sed//g'  > files.txt

while read i; do 
    python3 /scripts/get_ROH_demographics_bcftools_plink.py -i "$i".bed -o $i
	python3 /scripts/Transform_ROH_demographics.py -i ROH_demo_plots_"$i".txt -o ROH_demo_plots_"$i"
done 

cat *GP99_MAF_0.05_TV_1.6M_bcftools_plink_ROH*sum_transpose.txt | grep -v Samples | \
awk -F'\t' -v OFS='\t' '{print $1,"1.6M dataset",$2,$3}' > Combined_GP99_MAF_0.05_TV_1.6M_bcftools_plink_ROH_sum_transpose.txt

sed -i '1s/^/Sample\tDataset\tBins\tROH_SUM\n/' Combined_GP99_MAF_0.05_TV_1.6M_bcftools_plink_ROH_sum_transpose.txt

# combined dataset
# Filter bcftools ROH for only ROH >=4Mb and transform combined dataset into bed format
awk -v OFS='\t' '{ if ($6>=4000000) print $3,$4,$5,$2,$6}'  /published_paleogenomes_ROH/combined_glimpse2_GP99_MAF05_published_0-5X_geno0_new.refcheck.fix_10Q_200SNPs.txt \
 > combined_glimpse2_GP99_MAF05_published_0-5X_geno0_new.refcheck.fix_10Q_200SNPs_4mb.bed
 
# Transform Plink combined dataset to bed format
awk -v OFS='\t' '{print $4,$7,$8,$1,$9}' /published_paleogenomes_ROH/combined_glimpse2_GP99_MAF05_published_0-5X_geno0_500_100kb_1het_200SNP.hom \
 > combined_glimpse2_GP99_MAF05_published_0-5X_geno0_500_100kb_1het_200SNP.bed
 
# Combine and merge for each individual seperately
 
 cat combined_glimpse2_GP99_MAF05_published_0-5X_geno0_new.refcheck.fix_10Q_200SNPs_4mb.bed combined_glimpse2_GP99_MAF05_published_0-5X_geno0_500_100kb_1het_200SNP.bed \
 | grep -v CHR | awk -v OFS='\t' '{print $1,$2,$3,$4,$5}' | sort -k1,1n -k2,2n -k4,4 | awk '{print > $4".tmp.bed"}'
 
# Merge each ID group separately
for file in *.tmp.bed; do
    ID=$(basename "$file" .tmp.bed)
    bedtools merge -i "$file" -c 4,5 -o count,collapse
done > combined_glimpse2_GP99_MAF05_published_0-5X_geno0_bcftools_plink_ROH.bed

# Cleanup
rm *.tmp.bed

# calculate bin size + transform
python3 /scripts/get_ROH_demographics_bcftools_plink.py -i combined_glimpse2_GP99_MAF05_published_0-5X_geno0_bcftools_plink_ROH.bed \
-o combined_glimpse2_GP99_MAF05_published_0-5X_geno0_bcftools_plink_ROH
python3 /scripts/Transform_ROH_demographics.py -i ROH_demo_plots_combined_glimpse2_GP99_MAF05_published_0-5X_geno0_bcftools_plink_ROH.txt \
-o ROH_demo_plots_combined_glimpse2_GP99_MAF05_published_0-5X_geno0_bcftools_plink_ROH

awk -F'\t' -v OFS='\t' '{print $1,"Combined dataset",$2,$3}' ROH_demo_plots_combined_glimpse2_GP99_MAF05_published_0-5X_geno0_bcftools_plink_ROH_sum_transpose.txt \
 > Final_ROH_demo_plots_combined_glimpse2_GP99_MAF05_published_0-5X_geno0_bcftools_plink_ROH_sum_transpose.txt
 
 sed -i '1s/^/Sample\tDataset\tBins\tROH_SUM\n/' Final_ROH_demo_plots_combined_glimpse2_GP99_MAF05_published_0-5X_geno0_bcftools_plink_ROH_sum_transpose.txt
 
 # Combine for Figure 5
 
 cat Final_ROH_demo_plots_combined_glimpse2_GP99_MAF05_published_0-5X_geno0_bcftools_plink_ROH_sum_transpose.txt \
 Combined_GP99_MAF_0.05_TV_1.6M_bcftools_plink_ROH_sum_transpose.txt | grep -v Sample > Figure5_ROH_sum_profiles_plink_bcftools_comparison.txt
 
 sed -i '1s/^/Sample\tDataset\tBins\tROH_SUM\n/'  Figure5_ROH_sum_profiles_plink_bcftools_comparison.txt
