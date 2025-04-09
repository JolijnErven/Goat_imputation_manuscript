#!/bin/bash

# Script to add colour codes to ROH files for specific samples, 
# plot their local ROH profiles, and generate local heterozygosity plots.
# The script processes both PLINK and BCFtools-generated files, adds colour to them,
# and then uses R scripts to plot the ROH profiles and heterozygosity.

# Samples of interest
echo -e "Semnan3\nPotterne1\nBlagotin16\nBulak1\nGanjdareh22\nDirekli6\nKazbegi1" > samples.txt

# Get list of .hom files from PLINK for each sample and store in files_plink.txt
while read i; do ls *"$i"*.hom ; done < samples.txt > files_plink.txt

# Add colour to PLINK .hom files using a separate script
while read i; do /scripts/get_colour.sh $i ; done < files_plink.txt

# Repeat the process for BCFtools-generated files
while read i; do ls *"$i"*Q10_200SNPs.txt ; done < samples.txt > files_bcftools.txt

# Add colour to BCFtools ROH files using a separate script
while read i; do /scripts/get_colour_bcftools.sh $i ; done < files_bcftools.txt

# ROHan
while read i; do ls *"$i"*.max.hmmrohl.gz | sed 's/.gz//g' ; done < samples.txt > files_ROHan.txt

while read i; do gzip -d "$i".gz ; /scripts/get_colour_ROHan.sh $i ; done < files_ROHan.txt

# Plot local ROH profiles for each sample and chromosome
# Chromosomes of interest for each sample, listed in the same order as the samples
echo -e "15\n7\n9\n1\n1\n4\n24" > chrom_sam.txt

# Print the chromosome sizes associated with the chromosomes of interest
echo -e "81904557\n108433636\n91568626\n157403528\n157403528\n120734966\n62310066" > chrom_size.txt

# Create command for plotting by pasting together the sample, chromosome, and chromosome size information
paste samples chrom_sam.txt chrom_size.txt | awk '{ print "-s",$1,"-c",$2,"-e",$3}' > command.txt

# Run the local ROH profile plot R script for each sample, chromosome, and size
while read i; do /Software/R-4.2.0/bin/Rscript /plotting_scripts/Local_ROH_profile_plot.R $i ; done < command.txt

# Generate local heterozygosity plots for each sample and chromosome
# Extract relevant fields for heterozygosity plotting
cut -f 1,2,3,4 -d " " command.txt > het_command.txt

# Run the heterozygosity profile plot R script for each sample and chromosome
while read i; do /Software/R-4.2.0/bin/Rscript /plotting_scripts/Plot_heterozygosity_profile.R $i ; done < het_command.txt
