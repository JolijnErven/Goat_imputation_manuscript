#!/bin/bash

# This script calls variants from BAM files using SAMtools and BCFtools.
# It first checks if a reference site file exists; if not, it creates one.
# Then, it performs variant calling using mpileup/call and outputs a VCF file.
# Usage: run_variant_calling_filtering.sh <chromosome> <sample_id> <reference_vcf>

# Define directories, chr and samples
ref_dir="path_to_reference_panel_directory"  # Directory containing reference files
dir="path_to_bam_files"  # Directory containing BAM files
chr=$1  # Chromosome number (first argument)
sample=$2  # Sample ID (second argument, BAM file name)
reference_vcf=$3  # Reference VCF file name (third argument)
reference="/Reference_Genomes/goat/ARS1.fa"  # Reference genome

# Check if the reference site file and its index exist, if not, create them
if [ ! -f ${ref_dir}/chr${chr}_${reference_vcf}.sites.tsv.gz ] || [ ! -f ${ref_dir}/chr${chr}_${reference_vcf}.sites.tsv.gz.tbi ]; then
    echo "Creating reference site file for chr${chr}_${reference_vcf}.vcf.gz"
    
    # Extract the reference site information from the VCF and compress it
    /Software/bcftools-1.13/bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' \
    ${ref_dir}/chr${chr}_${reference_vcf}.vcf.gz | bgzip -c > ${ref_dir}/chr${chr}_${reference_vcf}.sites.tsv.gz
    
    # Index the site file using tabix
    tabix -s1 -b2 -e2 ${ref_dir}/chr${chr}_${reference_vcf}.sites.tsv.gz
    
    # Convert to BED format and compress with bgzip
    zcat ${ref_dir}/chr${chr}_${reference_vcf}.sites.tsv.gz | awk '{print $1"\t"$2-1"\t"$2}' | grep -v "#" | bgzip -c > ${ref_dir}/chr${chr}_${reference_vcf}.sites.bed
    tabix -p ${ref_dir}/chr${chr}_${reference_vcf}.sites.bed.gz
fi

# Perform variant calling using BCFtools mpileup and call
echo "Running BCFtools mpileup and variant calling for chr${chr} sample ${sample}..."

# Use samtools mpileup to create the pileup and BCFtools to call variants
/Software/samtools-1.13/samtools mpileup -f ${reference} -t SP,AD,INFO/AD,ADF,ADR,DP,INFO/DPR -B -q 30 -Q 30 -s -O -u \
-l ${ref_dir}/chr${chr}_${reference_vcf}.sites.bed.gz --threads 32 -r ${chr} ${dir}/${sample}.bam -Ou | \
/Software/bcftools-1.13/bcftools call -m -f GQ,GP -C alleles -T ${ref_dir}/chr${chr}_${reference_vcf}.sites.tsv.gz --threads 32 -Oz \
-o ${dir}/VCF/chr${chr}_${sample}_mpileup.vcf.gz

# Index the resulting VCF file
echo "Indexing the output VCF file..."
/Software/bcftools-1.13/bcftools index ${dir}/VCF/chr${chr}_${sample}_mpileup.vcf.gz

echo "Variant calling completed for chr${chr} sample ${sample}."

# Filtering steps for VCF
#Filter for depth and genotype quality
python /scripts/filter_hard_call_vcfs_GQ30-DP6-40.py ${dir}/VCF/chr${chr}_${sample}_mpileup.vcf.gz | bgzip -c > ${dir}/VCF/chr${chr}_${sample}_mpileup_DP6-40_GQ30.vcf.gz
