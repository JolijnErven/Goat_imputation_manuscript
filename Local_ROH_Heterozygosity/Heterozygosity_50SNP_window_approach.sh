#!/bin/bash

# Calculates heterozygosity in a 50-SNP window.
# Applies this both to individual and combined VCF files.

# Define samples of interest and save to file
echo -e 'Semnan3\nPotterne1\nBlagotin16\nBulak1\nGanjdareh22\nDirekli6\nKazbegi1' > samples.txt

# Get heterozygosity from individual vcf files (e.g. 1.6M and 1.2M)
while read i; do
    ls /published_paleogenomes_ROH/*"$i"*.vcf.gz | sed 's/.vcf.gz//g' | rev | cut -f 1 -d '/' | rev >> vcfs.txt
 done < samples.txt

while read vcf; do
    for c in {1..29}; do
        zcat /published_paleogenomes_ROH/"$vcf".vcf.gz | awk -v chr=${c} '{ if ($1==chr) print $0}' > chr"$c"_"$vcf".vcf ; awk -f divide_chr_50_SNPs.awk chr"$c"_"$vcf".vcf > hets_chr"$c"_"$vcf".vcf ; awk '{ if ($3>2) print $0}' hets_chr"$c"_"$vcf".vcf > figure_hets_chr"$c"_"$vcf".vcf
     done
 done < vcfs.txt

# Get heterozygosity from combined vcf files
while read s; do
    for c in {1..29}; do
	    bcftools view -s $s combined_glimpse2_GP99_MAF05_published_0-5X_geno0.refcheck.fix.vcf.gz -Ov -r $c -o chr"$c"_"$s"_combined_glimpse2_GP99_MAF05_published_0-5X_geno0.vcf
        awk -f divide_chr_50_SNPs.awk chr"$c"_"$s"_combined_glimpse2_GP99_MAF05_published_0-5X_geno0.vcf > hets_chr"$c"_"$s"_combined_glimpse2_GP99_MAF05_published_0-5X_geno0.vcf
		awk '{ if ($3>2) print $0}' hets_chr"$c"_"$s"_combined_glimpse2_GP99_MAF05_published_0-5X_geno0.vcf > figure_hets_chr"$c"_"$s"_combined_glimpse2_GP99_MAF05_published_0-5X_geno0.vcf
     done
 done < samples.txt
