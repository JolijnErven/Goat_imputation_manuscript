#!/bin/bash

# Pipeline: Refined IBD Analysis from Imputed Ancient Goat Genomes
# Steps:
#   1. Merge imputed VCFs
#   2. Run Refined IBD 3 times (to address phasing uncertainty)
#   3. Filter segments > 3 cM (then repeated for 4 cM)
#   4. Extract unique pairs of individuals with shared IBD
#   5. Extract IBD segments per pair
#   6. Merge overlapping IBD into BED format
#   7. Add genetic map (cM coordinates)
#   8. Classify relationship (related/unrelated)

# Step 3-8 are performed for a merge of 0.6cM and 4 cM

# STEP 1: Merge Imputed VCFs (GP99)
echo "Creating list of imputed VCFs..."
ls /Imputed/published/glimpse_Plink/*GP99.vcf.gz > merge_list.txt
ls /Imputed/glimpse_Plink/*0.5X*GP99.vcf.gz >> merge_list.txt

echo "Merging VCFs..."
/Software/bcftools-1.17/bcftools merge -l merge_list.txt -Oz -o Published_imputed_ancient_goat_genomes.vcf.gz
tabix Published_imputed_ancient_goat_genomes.vcf.gz

# STEP 2: Run Refined IBD (3 repetitions for potential phasing errors)
/scripts/run_refined_IBD_3_times.sh \
  Published_imputed_ancient_goat_genomes \
  Published_imputed_ancient_goat \
  _phased_glimpse2-filtering.vcf.gz \
  combined_above_3x_glimpse2_GP99_MAF_0.05_TV_no_low_conc_no_single_chr18_2mb_removed.pos.txt \
  100000

# STEP 3: Filter IBD segments > 3cM
echo "Filtering for segments longer than 3cM..."
awk '{ if ($9 > 3) print $0 }' \
  Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_0.6merge_refinedibd.ibd1.LOD3 \
  > Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_0.6merge_refinedibd.ibd1.LOD3.3cM

# STEP 4: Extract unique individual pairs
cut -f1,3 Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_0.6merge_refinedibd.ibd1.LOD3.3cM \
  | sort | uniq \
  > Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_0.6merge_refinedibd.ibd1.LOD3.3cM_pairs

# STEP 5: Generate commands to extract IBD segments per pair
awk -v dir="Published_imputed_ancient_goat_1" -v file="Published_imputed_ancient_goat_genomes" '{
  print "grep \"" "^" $1 "[[:space:]]\" " dir "/" file "_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_0.6merge_refinedibd.ibd1.LOD3.3cM" \
        " | grep \"" $2 "[[:space:]]\" > " dir "/" $1 "." $2 "." file "_MAF_0.05_no_missing_TV_0.6merge_refinedibd.ibd1.LOD3.3cM_duo"
}' Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_0.6merge_refinedibd.ibd1.LOD3.3cM_pairs \
> pull_pairs.lod3.sh

cat pull_pairs.lod3.sh | parallel -j15
wait

# STEP 6: Merge IBD segments into BED format
echo "Merging segments..."
for file in $(ls Published_imputed_ancient_goat_1/ | grep 0.6merge_| grep duo); do
  one=$(awk '{print $1}' Published_imputed_ancient_goat_1/$file | sort -u)
  two=$(awk '{print $3}' Published_imputed_ancient_goat_1/$file | sort -u)

  awk '{print $5"\t"$6"\t"$7"\t"$9}' Published_imputed_ancient_goat_1/$file \
    | sort -k1,1n -k2,2n \
    | mergeBed -c 4 -o max \
    | awk -v var="$one" -v yar="$two" '{print var, yar, $1, $2, $3, $4}' \
    > Published_imputed_ancient_goat_1/${one}-${two}.Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99__MAF_0.05_TV_0.6merge_LOD3.3cM.mergeBed
done

# STEP 7: Add cM coordinates from genetic map
echo "Adding genetic map (cM) coordinates..."
ls Published_imputed_ancient_goat_1/*0.6merge*mergeBed > bed_files.txt
while read b; do
  echo "/scripts/cm_add.sh $b"
done < bed_files.txt > parallel_cm_add.sh

cat parallel_cm_add.sh | parallel -j15
wait

# STEP 8: Compute relationship classification (related/unrelated)
echo "Classifying individual relationships..."
for file in $(ls Published_imputed_ancient_goat_1/ | grep 0.6merge | grep mergeBed.cm); do
  segs=$(wc -l Published_imputed_ancient_goat_1/$file | awk '{print $1}')
  sum=$(awk '{sum += $7 - $6} END {print sum}' Published_imputed_ancient_goat_1/$file)
  one=$(awk '{print $1}' Published_imputed_ancient_goat_1/$file | sort -u)
  two=$(awk '{print $2}' Published_imputed_ancient_goat_1/$file | sort -u)

  echo "$one $two $segs $sum" \
    | awk '{ if ($3 > 2 && $4 > 23) print $0, "rel"; else print $0, "un" }'
done > Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_TV_0.6merge.LOD3_cM.3cM_ind.txt


# repeat for 4cM
# STEP 3: Filter IBD segments > 3cM
echo "Filtering for segments longer than 3cM..."
awk '{ if ($9 > 3) print $0 }' \
  Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_4merge_refinedibd.ibd1.LOD3 \
  > Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_4merge_refinedibd.ibd1.LOD3.3cM

# STEP 4: Extract unique individual pairs
cut -f1,3 Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_4merge_refinedibd.ibd1.LOD3.3cM \
  | sort | uniq \
  > Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_4merge_refinedibd.ibd1.LOD3.3cM_pairs

# STEP 5: Generate commands to extract IBD segments per pair
awk -v dir="Published_imputed_ancient_goat_1" -v file="Published_imputed_ancient_goat_genomes" '{
  print "grep \"" "^" $1 "[[:space:]]\" " dir "/" file "_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_4merge_refinedibd.ibd1.LOD3.3cM" \
        " | grep \"" $2 "[[:space:]]\" > " dir "/" $1 "." $2 "." file "_MAF_0.05_no_missing_TV_4merge_refinedibd.ibd1.LOD3.3cM_duo"
}' Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_no_missing_TV_4merge_refinedibd.ibd1.LOD3.3cM_pairs \
> pull_pairs.lod3.sh

cat pull_pairs.lod3.sh | parallel -j15
wait

# STEP 6: Merge IBD segments into BED format
echo "Merging segments..."
for file in $(ls Published_imputed_ancient_goat_1/ | grep 4merge_| grep duo); do
  one=$(awk '{print $1}' Published_imputed_ancient_goat_1/$file | sort -u)
  two=$(awk '{print $3}' Published_imputed_ancient_goat_1/$file | sort -u)

  awk '{print $5"\t"$6"\t"$7"\t"$9}' Published_imputed_ancient_goat_1/$file \
    | sort -k1,1n -k2,2n \
    | mergeBed -c 4 -o max \
    | awk -v var="$one" -v yar="$two" '{print var, yar, $1, $2, $3, $4}' \
    > Published_imputed_ancient_goat_1/${one}-${two}.Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99__MAF_0.05_TV_4merge_LOD3.3cM.mergeBed
done

# STEP 7: Add cM coordinates from genetic map
echo "Adding genetic map (cM) coordinates..."
ls Published_imputed_ancient_goat_1/*4merge*mergeBed > bed_files.txt
while read b; do
  echo "/scripts/cm_add.sh $b"
done < bed_files.txt > parallel_cm_add.sh

cat parallel_cm_add.sh | parallel -j15
wait

# STEP 8: Compute relationship classification (related/unrelated)
echo "Classifying individual relationships..."
for file in $(ls Published_imputed_ancient_goat_1/ | grep 4merge | grep mergeBed.cm); do
  segs=$(wc -l Published_imputed_ancient_goat_1/$file | awk '{print $1}')
  sum=$(awk '{sum += $7 - $6} END {print sum}' Published_imputed_ancient_goat_1/$file)
  one=$(awk '{print $1}' Published_imputed_ancient_goat_1/$file | sort -u)
  two=$(awk '{print $2}' Published_imputed_ancient_goat_1/$file | sort -u)

  echo "$one $two $segs $sum" \
    | awk '{ if ($3 > 2 && $4 > 23) print $0, "rel"; else print $0, "un" }'
done > Published_imputed_ancient_goat_1/Published_imputed_ancient_goat_genomes_beagle_batch_impute_phase_GP99_MAF_0.05_TV_4merge.LOD3_cM.3cM_ind.txt
