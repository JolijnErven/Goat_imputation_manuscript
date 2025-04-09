#!/bin/bash

# ROH Exploration with PLINK for 0.5X Imputed Data
# -------------------------------------------------
# This script explores the minimum SNP threshold for a Runs of Homozygosity (ROH) analysis
# using PLINK on a dataset with 0.5X imputed samples and also explores the ideal number of sites.
# The analysis produced Figure S28 in the study.

# Set parameters
dir="/modern_reference_dataset/"
dataset="new_vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans"

# Step 1: Merge modern dataset with ancient imputed dataset (GP99 - 0.5X and high coverage imputed, and non-imputed high coverage samples)
plink --noweb --chr-set 29 \
      --bfile ${dir}${dataset} \
      --merge-list ancients.txt \
      --make-bed \
      --out ${dataset}_ancient_merged_GP99

# Step 2: Filter merged dataset for Minor Allele Frequency (MAF) >= 5% and no missingness
plink --bfile ${dataset}_ancient_merged \
      --chr-set 29 \
      --geno 0 \
      --make-bed \
      --extract ${dir}vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_MAF_0.05.pos.txt \
      --out ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0

# Step 3: Extract only imputed samples
grep Imp ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0.fam > keep.txt
mkdir -p SNP_exploration
plink --bfile ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0 \
      --chr-set 29 \
      --make-bed \
      --keep keep.txt \
      --out SNP_exploration/combined_${dataset}_0.5X_GP99_MAF_0.05_geno_0

plink --bfile ${dataset}_ancient_merged_GP99_MAF_0.05_geno_0 \
      --chr-set 29 \
      --make-bed \
      --keep keep.txt \
      --extract ${dir}vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_TV_MAF_0.05_pos.txt \
      --out SNP_exploration/combined_${dataset}_0.5X_GP99_MAF_0.05_geno_0_TV

# Step 4: Create individual bed files filtered for all SNP sites and transversions
# Get list of samples
ls /Imputed/glimpse_Plink/*0.5x*GP99.bed | sed 's/.bed//g' | rev | cut -f 1 -d '/' | rev > samples.txt
# Add high coverage imputed test samples 
ls /Imputed/glimpse_Plink/*HC*GP99.bed | sed 's/.bed//g'| rev | cut -f 1 -d '/' | rev >> samples.txt

# Generate individual bed files
while read i; do
  plink --bfile /Imputed/glimpse_Plink/"$i" \
        --chr-set 29 \
        --make-bed \
        --extract ${dir}vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_MAF_0.05.pos.txt \
        --out SNP_exploration/"$i"_MAF_0.05

  plink --bfile SNP_exploration/"$i"_MAF_0.05 \
        --chr-set 29 \
        --make-bed \
        --extract ${dir}vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_TV_MAF_0.05_pos.txt \
        --out SNP_exploration/"$i"_MAF_0.05_TV

done < samples.txt

# Step 6: Perform SNP and site exploration for each bed file in SNP_exploration directory
echo -e "1000000\n750000\n500000" > SNP_exploration/downsample.txt
cd SNP_exploration/
ls *bed | sed 's/.bed//g' > samples.txt

while read i; do 
  cut -f 2  "$i".bim > "$i"_pos.txt
  
  for s in 50 100 150 200 250 300; do
    plink --bfile "$i" \
          --chr-set 29 \
          --homozyg \
          --homozyg-density 50 \
          --homozyg-gap 100 \
          --homozyg-kb 500 \
          --homozyg-snp $s \
          --homozyg-window-het 1 \
          --homozyg-window-snp 50 \
          --homozyg-window-threshold 0.05 \
          --memory 40000 \
          --threads 4 \
          --out "$i"_Full_500kb_100kb_1het_"$s"SNPs
  done
  
  while read d; do 
    shuf -n $d "$i"_pos.txt | sort -t: -k1,1 -k2,2 -n > "$i"_"$d".txt
    plink --bfile "$i" \
          --extract "$i"_"$d".txt \
          --make-bed \
          --chr-set 29 \
          --out "$i"_"$d" \
          --threads 4 \
          --memory 40000
    
    for s in 50 100 150 200 250 300; do
      plink --bfile "$i"_"$d" \
            --chr-set 29 \
            --homozyg \
            --homozyg-density 50 \
            --homozyg-gap 100 \
            --homozyg-kb 500 \
            --homozyg-snp $s \
            --homozyg-window-het 1 \
            --homozyg-window-snp 50 \
            --homozyg-window-threshold 0.05 \
            --memory 40000 \
            --threads 4 \
            --out "$i"_"$d"_500kb_100kb_1het_"$s"SNPs
    done
  done < downsample.txt

done < samples.txt

# Step 6: Get individual high coverage non imputed for comparison (not downsampled)
ls /high_coverage/*bed | sed 's/.bed//g' | rev | cut -f 1 -d '/' | rev > HC_samples.txt

while read i; do
  plink --bfile /high_coverage/"$i" \
        --chr-set 29 \
        --make-bed \
        --extract ${dir}vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_MAF_0.05.pos.txt \
        --out "$i"_HC_MAF_0.05

  plink --bfile "$i"_HC_MAF_0.05 \
        --chr-set 29 \
        --make-bed \
        --extract ${dir}vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_TV_MAF_0.05_pos.txt \
        --out "$i"_HC_MAF_0.05_TV

done < HC_samples.txt

# Step 6.1 Get ROH from high coverage using default parameters

ls *HC_MAF*bed | sed 's/.bed//g' > HC_samples.txt

while read i; do
plink --bfile "$i" \
    --chr-set 29 \
    --homozyg \
    --homozyg-density 50 \
    --homozyg-gap 100 \
    --homozyg-kb 500 \
    --homozyg-snp 50 \
    --homozyg-window-het 1 \
    --homozyg-window-snp 50 \
    --homozyg-window-threshold 0.05 \
    --memory 40000 \
    --threads 4 \
    --out "$i"_500kb_100kb_1het_50SNPs
; done < HC_samples.txt
	
# Step 7: Transform .hom file into ROH bin sizes and ease for plotting

ls *.hom | sed 's/.hom//g' > hom.samples.txt

while read i; do 
    python3 /scripts/get_ROH_demographics.py -i "$i".hom -o $i
	
	python3 /scripts/Transform_ROH_demographics.py -i ROH_demo_plots_"$i".txt -o ROH_demo_plots_"$i"

; done < hom.samples.txt

# Step 8: Add SNP, downsample and dataset information to output files
# Step 8.1 Add SNP and downsample information
# Add Full to downsample for going over iterations
echo "Full" >> downsample.txt

while read i; do 
    while read d; do 
        for s in 50 100 150 200 250 300; do
            awk -F'\t' -v OFS='\t' -v snp=$s -v down=$d '{print $1_snp"SNPs", down, $2, $3}' \
            ROH_demo_plots_"$i"_"$d"_500kb_100kb_1het_"$s"SNPs_sum_transpose.txt \
            > ROH_demo_plots_"$i"_"$d"_500kb_100kb_1het_"$s"SNPs_sum_transpose_info.txt
        done 
    done < downsample.txt
done < sample.txt

# Step 8.2: Add dataset information
ls *combined_*sum_transpose_info.txt | grep -v _TV_ > combined_dataset_all_sites.txt
while read i; do
    awk -F'\t' -v OFS='\t' '{print $1,$2,"Combined dataset","All sites",$3,$4}' $i > Final_"$i"
; done < combined_dataset_all_sites.txt

ls *combined_*_TV_*sum_transpose_info.txt  > combined_dataset_TV.txt
while read i; do
    awk -F'\t' -v OFS='\t' '{print $1,$2,"Combined dataset","Transversions",$3,$4}' $i > Final_"$i"
; done < combined_dataset_all_sites.txt

ls *sum_transpose_info.txt | grep -v _TV_ | grep -v combined > ind_dataset_all_sites.txt
while read i; do
    awk -F'\t' -v OFS='\t' '{print $1,$2,"Individual dataset","All sites",$3,$4}' $i > Final_"$i"
; done < combined_dataset_all_sites.txt

ls *_TV_*sum_transpose_info.txt | grep -v combined > ind_dataset_TV.txt
while read i; do
    awk -F'\t' -v OFS='\t' '{print $1,$2,"individual dataset","Transversions",$3,$4}' $i > Final_"$i"
; done < combined_dataset_all_sites.txt

# Add downsample information to non-downsampled high coverage for ease of visualization 
ls *HC_MAF_0.05*sum_transpose.txt | sed 's/.txt//g' > high_coverage_samples.txt

while read i; do 
    while read d; do 
         awk -F'\t' -v OFS='\t' -v down=$d '{print $1, down, $2, $3}' \
            "$i".txt  > "$i"_info.txt
	; done < downsample.txt
; done < high_coverage_samples.txt

grep MAF_0.05_TV high_coverage_samples.txt > high_coverage_samples_TV.txt

grep -v TV high_coverage_samples.txt > high_coverage_samples_all.txt

while read i; do
    awk -F'\t' -v OFS='\t' '{print $1, $2, "Individual dataset", "All sites", $3, $4}' "$i"_info.txt > Final_"$i"_info_ind.txt
    awk -F'\t' -v OFS='\t' '{print $1, $2, "Combined dataset", "All sites", $3, $4}' "$i"_info.txt > Final_"$i"_info_combined.txt
done < high_coverage_samples_all.txt
	
while read i; do
    awk -F'\t' -v OFS='\t' '{print $1, $2, "Individual dataset", "Transversions", $3, $4}' "$i"_info.txt > Final_"$i"_info_ind.txt
    awk -F'\t' -v OFS='\t' '{print $1, $2, "Combined dataset", "Transversions", $3, $4}' "$i"_info.txt > Final_"$i"_info_combined.txt
done < high_coverage_samples_TV.txt

# Step 9: Cat all seperate ROH sum statistic files
cat Final*sum_transpose_info*txt | grep -v Samples > Combined_ROH_sum_statistics_SNP_site_exploration_ready_for_plot.txt

# Add proper header
sed -i '1s/^/Sample\tDownsampled_sites\tDataset\tSites\tBins\tROH_SUM\n/' Combined_ROH_sum_statistics_SNP_site_exploration_ready_for_plot.txt

# Remove individual statistic files
rm ROH_demo*sum_transpose*txt