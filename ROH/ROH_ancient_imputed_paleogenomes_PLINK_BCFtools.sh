#!/bin/bash

# ROH Exploration with PLINK for >=0.5X Imputed previously published paleogenomes
# Tested individuals together with no missingness and individual-based calculations on 1.2 and 1.6 million sites
# The analysis produced Figure 5 and S32 in the study.

# Merge all imputed previously published >=0.5X paleogenomes, filter for MAF 5% and no missingness
mkdir -p published_paleogenomes_ROH
cd published_paleogenomes_ROH
ls /Imputed/published/glimpse_Plink/*GP99*bed | sed 's/.bed//g' > samples.txt
# Get first sample
first=$(head -1 samples.txt)
# Remove first sample from file
sed -i '1d' samples.txt

# Merge with plink
plink --bfile ${first} --chr-set 29 --geno 0 --merge-list samples.txt --make-bed
   --extract /modern_reference_dataset/vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_MAF_0.05.pos.txt
   --out published_paleogenomes_MAF_0.05_no_missingness
   
# ROH calculation with PLINK with a minimum of 200 SNPs
 plink --bfile published_paleogenomes_MAF_0.05_no_missingness \
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
    --out published_paleogenomes_MAF_0.05_no_missingness_500kb_100kb_1het_200SNPs
	
# Get ROH Bin statistics and transpose
python3 /scripts/get_ROH_demographics.py -i published_paleogenomes_MAF_0.05_no_missingness_500kb_100kb_1het_200SNPs.hom \
-o published_paleogenomes_MAF_0.05_no_missingness_500kb_100kb_1het_200SNPs
python3 /scripts/Transform_ROH_demographics.py -i ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness_500kb_100kb_1het_200SNPs.txt \
-o ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness_500kb_100kb_1het_200SNPs

# Add details to ROH sum transpose file
awk -F'\t' -v OFS='\t' '{print "$1","Combined dataset" $2, $3}' ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness_500kb_100kb_1het_200SNPs_sum_transpose.txt \
> ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness_500kb_100kb_1het_200SNPs_sum_transpose_info.txt

# Create individual 1.2 and 1.6m bed files and calculate ROH with a minimum of 200 SNPs
ls /Imputed/published/glimpse_Plink/*GP99*bed | rev | cut -f 1 -d '/'| rev | sed 's/.bed//g' > samples.txt

while read i; do 
plink --bfile /Imputed/published/glimpse_Plink/$i \
        --chr-set 29 \
        --make-bed \
        --extract /modern_reference_dataset/vargoats_snps_filtered_1270_beagle5-3_combined_updatedIDs_hans_TV_MAF_0.05_pos.txt \
        --out "$i"_MAF_0.05_TV
		
cut -f 2 "$i"_MAF_0.05_TV.bim > "$i"_MAF_0.05_TV_pos.txt
shuf -n 1200000 "$i"_MAF_0.05_TV_pos.txt | sort -t: -k1,1 -k 2,2 -n > "$i"_MAF_0.05_TV_1.2M.txt

plink --bfile "$i"_MAF_0.05_TV --extract "$i"_MAF_0.05_TV_1.2M.txt --make-bed  --chr-set 29 --out "$i"_MAF_0.05_TV_1.2M --threads 12 --memory 40000

plink --bfile "$i_MAF_0.05_TV"_1.2M  --chr-set 29 --homozyg --homozyg-density 50 --homozyg-gap 100 --homozyg-kb 500 --homozyg-snp 200 --homozyg-window-het 1 \
--homozyg-window-snp 50 --homozyg-window-threshold 0.05 --memory 40000 --threads 2 --out "$i"_MAF_0.05_TV_1.2M_500kb_100kb_200SNPs_1het


shuf -n 1600000 "$i"_MAF_0.05_TV_pos.txt | sort -t: -k1,1 -k 2,2 -n > "$i"_MAF_0.05_TV_1.6M.txt

plink --bfile "$i"_MAF_0.05_TV --extract "$i"_MAF_0.05_TV_1.6M.txt --make-bed  --chr-set 29 --out "$i"_MAF_0.05_TV_1.6M --threads 12 --memory 40000

plink --bfile "$i"_MAF_0.05_TV_1.6M  --chr-set 29 --homozyg --homozyg-density 50 --homozyg-gap 100 --homozyg-kb 500 --homozyg-snp 200 \
--homozyg-window-het 1 --homozyg-window-snp 50 --homozyg-window-threshold 0.05 --memory 40000 --threads 2 --out "$i"_MAF_0.05_TV_1.6M_500kb_100kb_200SNPs_1het

; done < samples.txt

# Get ROH Bin statistics and transpose
ls *hom | grep -E "1.2M|1.6M" | sed 's/.hom//g' > files.txt

while read i; do 
    python3 /scripts/get_ROH_demographics.py -i "$i".hom -o $i
	
	python3 /scripts/Transform_ROH_demographics.py -i ROH_demo_plots_"$i".txt -o ROH_demo_plots_"$i"

; done < files.txt

# Cat statistic files and add details to ROH sum transpose files
cat ROH_demo_plots*1.2M*sum_transpose.txt | grep -v variable > Combined_imputed_ancients_MAF_0.05_TV_1.2M_500kb_100kb_200SNPs_1het_sum_transpose.txt
cat ROH_demo_plots*1.6M*sum_transpose.txt | grep -v variable > Combined_imputed_ancients_MAF_0.05_TV_1.6M_500kb_100kb_200SNPs_1het_sum_transpose.txt

awk -F'\t' -v OFS='\t' '{print "$1","Individual 1.2M Transversion" $2, $3}' Combined_imputed_ancients_MAF_0.05_TV_1.2M_500kb_100kb_200SNPs_1het_sum_transpose.txt \
> Combined_imputed_ancients_MAF_0.05_TV_1.2M_500kb_100kb_200SNPs_1het_sum_transpose_info.txt
awk -F'\t' -v OFS='\t' '{print "$1","Individual 1.6M Transversion" $2, $3}' Combined_imputed_ancients_MAF_0.05_TV_1.6M_500kb_100kb_200SNPs_1het_sum_transpose.txt \
> Combined_imputed_ancients_MAF_0.05_TV_1.6M_500kb_100kb_200SNPs_1het_sum_transpose_info.txt


# BCFtools section
# Calculate ROH on combined dataset and 1.6M dataset used previously for plink
# First get VCF file
plink --bfile published_paleogenomes_MAF_0.05_no_missingness --chr-set 29  --recode vcf
--out published_paleogenomes_MAF_0.05_no_missingness --double-id

bgzip published_paleogenomes_MAF_0.05_no_missingness.vcf
tabix published_paleogenomes_MAF_0.05_no_missingness.vcf.gz

# Fix alleles in vcf (Due to plink switching alleles)
/Software/bcftools-1.17/bcftools norm -f Reference_Genomes/goat/ARS1.fa \
    --check-ref ws --do-not-normalize -Oz -o published_paleogenomes_MAF_0.05_no_missingness.refcheck.vcf.gz \
	published_paleogenomes_MAF_0.05_no_missingness.vcf.gz

zcat published_paleogenomes_MAF_0.05_no_missingness.refcheck.vcf.gz | grep "#" > published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix.vcf
zcat published_paleogenomes_MAF_0.05_no_missingness.refcheck.vcf.gz | awk '$4 != $5' | grep -v "#" >> published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix.vcf
bgzip published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix.vcf
tabix published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix.vcf.gz

# Calculate ROHs with bcftools
/Software/bcftools-1.17/bcftools roh -G 30 --AF-dflt 0.4 published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix.vcf.gz published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix.txt
# Filter ROH for minmum 200SNPs, 10 quality and a length of 500kb 
grep RG published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix.txt | awk '{if ($8>10) print $0}' |  awk '{if ($7>=200) print $0}' \
 | awk '{if ($6>=500000) print $0}' > published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q.txt

# Transform .txt file into ROH bin sizes and ease for plotting (Figure associated with output of this step is on OSF)
python3 /scripts/get_ROH_demographics_BCFtools.py -i published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q.txt \
-o published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q
python3 /scripts/Transform_ROH_demographics.py -i ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q.txt \
-o ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q

awk -F'\t' -v OFS='\t' '{print "$1","BCFtools combined" $2, $3}' ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q_sum_transpose.txt \
 | grep -v Sample > ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q_sum_transpose_info.txt

# Calculate ROH for 1.6M sites dataset on the same sites as plink
while read i; do 
# Get positions and extract them from imputed VCF
cut -f 1,4 "$i"_MAF_0.05_TV_1.6M.bim > "$i"_MAF_0.05_TV_1.6M.pos.txt
bcftools view -R "$i"_MAF_0.05_TV_1.6M.pos.txt /Imputed/published/GLIMPSE_ligate/"$i".vcf.gz -Oz -o "$i"_MAF_0.05_TV_1.6M.vcf.gz --threads 20
tabix "$i"_MAF_0.05_TV_1.6M.vcf.gz
# Calculate ROH
/raid_md0/jolijn/Software/bcftools-1.17/bcftools roh -G 30 --AF-dflt 0.4 "$i"_MAF_0.05_TV_1.6M.vcf.gz -o "$i"_MAF_0.05_TV_1.6M.txt
# Filter ROH for minmum 200SNPs, 10 quality and a length of 500kb 
grep RG "$i"_MAF_0.05_TV_1.6M.txt | awk '{if ($8>10) print $0}' |  awk '{if ($7>=200) print $0}' \
 | awk '{if ($6>=500000) print $0}' > "$i"_MAF_0.05_TV_1.6M_200SNPs_10Q.txt
# Transform .hom file into ROH bin sizes and ease for plotting 
    python3 python3 /scripts/get_ROH_demographics_BCFtools.py -i "$i"_MAF_0.05_TV_1.6M_200SNPs_10Q.txt -o "$i"_MAF_0.05_TV_1.6M_200SNPs_10Q
	python3 /scripts/Transform_ROH_demographics.py -i ROH_demo_plots_"$i"_MAF_0.05_TV_1.6M_200SNPs_10Q.txt -o ROH_demo_plots_"$i"_MAF_0.05_TV_1.6M_200SNPs_10Q
; done < samples.txt

# Cat statistic files and add details to ROH sum transpose files
cat ROH_demo_plots_*_MAF_0.05_TV_1.6M_200SNPs_10Q_sum_transpose.txt | grep -v variable > Combined_imputed_ancients_MAF_0.05_TV_1.6M_00SNPs_10Q_sum_transpose.txt
awk -F'\t' -v OFS='\t' '{print "$1","BCFtools Individual 1.6M Transversion" $2, $3}' Combined_imputed_ancients_MAF_0.05_TV_1.6M_00SNPs_10Q_sum_transpose.txt \
 > Combined_imputed_ancients_MAF_0.05_TV_1.6M_200SNPs_10Q_sum_transpose_info.txt
 
# Add header
sed -i '1s/^/Sample\tDataset\tBins\tROH_SUM\n/' ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness_500kb_100kb_1het_200SNPs_sum_transpose_info.txt
sed -i '1s/^/Sample\tDataset\tBins\tROH_SUM\n/' Combined_imputed_ancients_MAF_0.05_TV_1.2M_500kb_100kb_200SNPs_1het_sum_transpose_info.txt
sed -i '1s/^/Sample\tDataset\tBins\tROH_SUM\n/' Combined_imputed_ancients_MAF_0.05_TV_1.6M_500kb_100kb_200SNPs_1het_sum_transpose_info.txt
sed -i '1s/^/Sample\tDataset\tBins\tROH_SUM\n/' Combined_imputed_ancients_MAF_0.05_TV_1.6M_200SNPs_10Q_sum_transpose_info.txt
sed -i '1s/^/Sample\tDataset\tBins\tROH_SUM\n/' ROH_demo_plots_published_paleogenomes_MAF_0.05_no_missingness.refcheck.fix_200SNPs_10Q_sum_transpose_info.txt

# Remove individual files
rm ROH_demo_plots*1.6M*sum_transpose*txt
rm ROH_demo_plots*1.2M*sum_transpose*txt