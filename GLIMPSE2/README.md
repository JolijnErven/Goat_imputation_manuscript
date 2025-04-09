# GLIMPSE2 Imputation Pipeline

This directory contains scripts to perform genotype imputation using **GLIMPSE2** on downsampled BAM files, followed by filtering on genotype probabilities (GP) and conversion to PLINK BED format.

The imputation process relies on BAMs created in the `downsampling_BAM` directory.

---

## Step 1: Genotype Imputation with GLIMPSE2

GLIMPSE2 is used for imputation of genotypes from low-coverage sequencing data. The script `Run_GLIMPSE2.sh` orchestrates the imputation process.


### Reference
This step follows the official GLIMPSE2 pipeline:
- [GLIMPSE2 Documentation](https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/)

Ensure you have reference panels, a genetic map (optional), and the required GLIMPSE2 tools installed and accessible in your environment.

---

## Step 2: Filter by GP, Merge VCFs, and Convert to BED

After imputation, VCFs are filtered based on genotype probabilities (GP). The filtered files are then concatenated and transformed into PLINK BED format using the `Filter_GPs_merge_transform_to_bed.sh` script.


## Notes
- Input BAMs must be generated beforehand in the `downsampling_BAM` directory.
- Modify file paths in the scripts as needed to match your directory structure.
