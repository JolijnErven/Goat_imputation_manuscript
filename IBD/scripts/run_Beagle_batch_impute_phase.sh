#!/bin/bash

# Input Parameters:
# $1 = Input VCF prefix (expects ${1}.vcf.gz)
# $2 = Chromosome number (1–29)
# $3 = Modern reference panel
# $4 = Output directory
# $5 = Effective population size (ne)

# Beagle phasing and imputation
/Software/java8/bin/java -Xmx250g -jar /Software/Beagle/beagle.28Jun21.220.jar \
    gt=${1}.vcf.gz ref=/modern_reference_dataset/vargoats_snps_filtered_1372_beagle5-3_chr${2}_${3} \
    chrom=${2} impute=true window=40 overlap=4 gp=true seed=$RANDOM \
    out=${4}/chr${2}_${1}_beagle_batch_impute_phase \
    map=/raid_md0/goat/vargoats/bertrand_servin_map_breed_sex/Recombination_maps_1Mb_global_population_breed_sex-averaged_chr${2}.gmap \
    nthreads=40 ne=${5}
