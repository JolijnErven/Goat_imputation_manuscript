#!/bin/python

# Updated version of Rui Martiniano's script. Updated by Kevin Daly
# Script to filter phased genotypes with genotype probability lower than the threshold.
# This script also replaces triallelic genotypes with "./." and removes SNP counters.
# It can process both VCF and gzipped VCF files.
# Usage: python filter_vcf_GP.py -i input.vcf -o output.vcf -t 0.90

from __future__ import print_function
import argparse
import sys
import time
import gzip

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Filters phased genotypes with genotype probability lower than the threshold provided. It also replaces triallelic genotypes with './.'")

# Example usage:
# python filter_vcf_GP.py -i chr1.SAMPLE.vcf -o chr1.SAMPLE_GP90.vcf -t 0.90

# Define command-line arguments
parser.add_argument('-i', action="store", dest="vcf_input", type=str, help="Input VCF file (can be .vcf or .vcf.gz)")
parser.add_argument('-o', action="store", dest="vcf_output", type=str, help="Output VCF file")
parser.add_argument('-t', action="store", dest="genotype_probability_threshold", type=float, help="Threshold for genotype probability filtering")


try:
    options = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

# Get the input/output files and threshold from arguments
vcf_input = options.vcf_input
vcf_output = options.vcf_output
threshold = options.genotype_probability_threshold

# Open the output file for writing
outfile = open(vcf_output, 'w')

def genotype_prob_parser(snp, threshold, GP_index):
    """
    Filters genotypes with genotype probability lower than the threshold.
    Also replaces triallelic genotypes with './.'.
    """
    changed_snps = 0
    snps2 = 0
    empty_individual_snp = "./."
    new_snp = []
    
    for individual_snp in snp:
        genotype = individual_snp.split(":")[0]
        GP = individual_snp.split(":")[GP_index]
        
        # Filter based on genotype type and GP value
        if genotype.count("0") == 2:  # Homozygous reference
            if float(GP.split(",")[0]) >= threshold:
                new_snp.append(individual_snp)
            else:
                new_snp.append(empty_individual_snp)
                changed_snps += 1
        elif genotype.count("0") == 1:  # Heterozygous
            if float(GP.split(",")[1]) >= threshold:
                new_snp.append(individual_snp)
            else:
                new_snp.append(empty_individual_snp)
                changed_snps += 1
        elif genotype.count("1") == 2:  # Homozygous alternate
            if float(GP.split(",")[2]) >= threshold:
                new_snp.append(individual_snp)
            else:
                new_snp.append(empty_individual_snp)
                changed_snps += 1
        else:  # Triallelic genotypes
            new_snp.append(empty_individual_snp)
    
    return new_snp

# Process the VCF input file
counter = 0
header = []

# If the VCF input file is plain text (.vcf)
if vcf_input.endswith("vcf"):
    with open(vcf_input, 'r') as f:
        for line in f:
            # Process header lines
            if line.startswith('#'):
                headline = line.strip('\n').split("\t")
                header.append(headline)
                outfile.write('\t'.join(headline) + '\n')
            else:
                counter += 1
                # Process SNP lines with the genotype_prob_parser function
                info = line.strip('\n').split("\t")[0:9]
                format = line.strip('\n').split("\t")[8]
                GP_index = format.split(":").index("GP")
                snp = line.strip('\n').split("\t")[9:]
                outfile.write('\t'.join(info) + '\t' + '\t'.join(genotype_prob_parser(snp, threshold, GP_index)) + '\n')
    
    outfile.close()

# If the VCF input file is gzipped (.vcf.gz)
elif vcf_input.endswith("gz"):
    with gzip.open(vcf_input, 'rb') as f:
        for line in f:
            # Process header lines
            if line.startswith('#'):
                headline = line.strip('\n').split("\t")
                header.append(headline)
                outfile.write('\t'.join(headline) + '\n')
            else:
                counter += 1
                # Process SNP lines with the genotype_prob_parser function
                info = line.strip('\n').split("\t")[0:9]
                format = line.strip('\n').split("\t")[8]
                GP_index = format.split(":").index("GP")
                snp = line.strip('\n').split("\t")[9:]
                outfile.write('\t'.join(info) + '\t' + '\t'.join(genotype_prob_parser(snp, threshold, GP_index)) + '\n')
    
    outfile.close()

else:
    print("Please provide a .vcf or .vcf.gz file as input")
