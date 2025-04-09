#!/bin/bash

# Input variables:
# $1 = Prefix of input VCF
# $2 = Output prefix for directory
# $3 = Modern reference panel (VCF)
# $4 = Positions file to keep
# $5 = Effective population size (Ne) for phasing with beagle

# Step 1: Create directories for 3 replicates
for d in {1..3}; do
    mkdir -p ${2}_${d}
done

# Step 2: Phase dataset (Beagle) for all chromosomes in each replicate
for d in {1..3}; do
    for c in {1..29}; do
        /script/run_Beagle_batch_impute_phase.sh $1 $c $3 ${2}_${d} $5 
    done
done

# Step 3: Set genotype probabilities (GP) < 0.99 to missing
for d in {1..3}; do
    for c in {1..29}; do
        input_vcf="${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase.vcf.gz"
        output_vcf="${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase_GP99.vcf.gz"
        
        tabix "$input_vcf"
        /Software/bcftools-1.17/bcftools +setGT -Oz -o "$output_vcf" -- "$input_vcf" \
            -t q -i 'SMPL_MAX(FORMAT/GP) < 0.99' -n './.'
        tabix "$output_vcf"
    done
done

# Step 4: Filter sites with missing genotypes + keep only positions in $4
for d in {1..3}; do
    for c in {1..29}; do
        filtered_vcf="${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase_GP99_MAF_0.05.vcf.gz"
        out_vcf="${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV.vcf.gz"
        
        tabix "$filtered_vcf"
        vcftools --gzvcf "$filtered_vcf" --positions $4 --max-missing 1 \
            --recode --recode-INFO-all --stdout | bgzip -c > "$out_vcf"
    done
done

# Step 5: Tabix index all final VCFs before IBD
for d in {1..3}; do
    for c in {1..29}; do
        tabix "${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase_GP99_MAF_0.025_no_missing_TV.vcf.gz"
    done
done

# Step 6: Prepare Refined IBD commands and run in parallel
for d in {1..3}; do
    for c in {1..29}; do
        echo "java -Xss50m -jar /raid_md0/Software/refined-ibd.17Jan20.102.jar \
gt=${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse.vcf.gz \
out=${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse_refinedIBD \
nthreads=4 chrom=${c} \
map=/raid_md0/goat/vargoats/bertrand_servin_map_breed_sex/Recombination_maps_1Mb_global_population_breed_sex-averaged_chr${c}.gmap" \
        >> run_ibd.sh
    done
done

cat run_ibd.sh | parallel -j10
wait


# Step 7: Merge Refined IBD segments with two parameter sets
for d in {1..3}; do
    for c in {1..29}; do
        input_ibd="${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse_refinedIBD.ibd.gz"
        vcf_file="${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse.vcf.gz"
        map_file="/vargoats/bertrand_servin_map_breed_sex/Recombination_maps_1Mb_global_population_breed_sex-averaged_chr${c}.gmap"

        # Merge with 0.6cM gap, 1 discordant homozygote
        zcat "$input_ibd" | java -jar /raid_md0/jolijn/Software/merge-ibd-segments.17Jan20.102.jar \
            "$vcf_file" "$map_file" 0.6 1 \
            > "${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse.refinedIBD.ibd.gap_0.6_merge"

        # Merge with 4cM gap, 40 discordant homozygotes
        zcat "$input_ibd" | java -jar /raid_md0/jolijn/Software/merge-ibd-segments.17Jan20.102.jar \
            "$vcf_file" "$map_file" 4 40 \
            > "${2}_${d}/chr${c}_${1}_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse.refinedIBD.ibd.gap_4_merge"
    done
done


# Step 8: Combine all merged IBD segments (across replicates) and filter for LOD > 3
# Do this for 0.6 cM gap merge
cat ${2}_1/*${1}*_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse*gap_0.6_merge \
 ${2}_2/*${1}*_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse*gap_0.6_merge \
 ${2}_3/*${1}*_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse*gap_0.6_merge \
 | awk '$8 > 3' > ${2}_1/Combined_${1}_glimpse_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_0.6merge_refinedibd.ibd1.LOD3
# and 4 cM gap merge
cat ${2}_1/*${1}*beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse*gap_4_merge \
 ${2}_2/*${1}*beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse*gap_4_merge \
 ${2}_3/*${1}*beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_glimpse*gap_4_merge \
 | awk '$8 > 3' > ${2}_1/Combined_${1}_glimpse_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_4merge_refinedibd.ibd1.LOD3
