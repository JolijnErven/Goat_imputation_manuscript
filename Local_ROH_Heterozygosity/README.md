# Local ROH and Heterozygosity Profiles Pipeline

This directory contains scripts for generating local Runs of Homozygosity (ROH) and heterozygosity profiles across the genome. These profiles provide insight into fine-scale patterns of genetic diversity.

---

## Overview
The pipeline consists of two main steps:

### Step 1: Calculate Heterozygosity in 50-SNP Windows

This step computes heterozygosity using non-overlapping windows of 50 SNPs.

#### Script
```
Heterozygosity_50SNP_window_approach.sh
```

- Calculates the number of heterozygous calls in each 50-SNP window
- Requires `.vcf.gz` files generated from previous ROH analysis (see ROH pipeline)

### Step 2: Generate and Plot Local ROH and Heterozygosity Profiles

This step enhances ROH profile files with additional annotations and generates visualizations of both heterozygosity and ROH signals.

#### Script
```
Get_local_ROH_profiles.sh
```

- Requires ROH files generated from previous ROH analysis (see ROH pipeline)
- Plots local ROH profiles and local heterozygosity

---

## Input Requirements
- ROH files from prior ROH calling (e.g., PLINK (.hom) and bcftools (10Q_200SNPs.txt))
- VCF or genotype files compatible with heterozygosity script (Output from Heterozygosity_50SNP_window_approach.sh)

## Notes
- Ensure all required tools (e.g., `plink`, `bedtools`, R, etc.) are installed
- For large genomes or many samples, use job scheduling or parallelization as needed


