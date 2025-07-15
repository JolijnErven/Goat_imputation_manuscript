#!/bin/bash

# Filter imputed VCFs for GP, merge, transform to BED and remove low concordance regions


# Step 1: Filter imputed VCFs Based on GP
cd /Imputed/glimpse_Ligated/
echo -e '0.70\n0.80\n0.90\n0.95\n0.99' > GPs.txt
ls /*vcf.gz | sed 's/.vcf.gz//g' > vcfs.txt

while read i; do
    while read g; do
	    GP = $(echo $g | cut -f 2 -d '.')
	    python2 /scripts/filter_vcf_GP_v2.py "$i".vcf.gz -o "$i"_GP${GP}.vcf -t $g 
		bgzip "$i"_GP${GP}.vcf 
		tabix "$i"_GP${GP}.vcf.gz
	; done < GPs.txt
; done < vcfs.txt

# Step 2: Concatenate filtered VCFs
ls chr9_*vcf.gz | cut -f 2-100 -d "_" > concat_vcfs.txt

while read suffix; do
    files=""
    for chr in {1..29}; do
        files+="Imputed.chr${chr}_${suffix} "
    done
    /Software/bcftools-1.12/bcftools concat $files -Oz -o Merged_${suffix}.vcf.gz --threads 20
done < concat_vcfs.txt

# Step 3: Convert to PLINK BED format and exclude low concordance and FPR positions
ls Merged*vcf.gz | sed 's/.vcf.gz//g' > to_bed.txt

while read i; do
    plink --vcf "$i".vcf.gz --chr-set 29 --make-bed --exclude /validation/local/Positions_to_remove_low_concordance_FPR.pos.txt \
	--out ../glimpse_Plink/"$i"
; done 	< to_bed.txt

