# Sliding Window Concordance and FPR Analysis Pipeline

This directory contains scripts to calculate non reference concordance and heterozygous false positive rates (FPR) in sliding windows across the genome. It also generates Manhattan plots to visualize these metrics and identify problematic regions for downstream analysis.

---

## Step 1: Identify Matching and Discordant Genotypes

To compute local concordance and FPR, we first identify matching and discordant genotypes using the following R script:

### Script
```
Concordance_FPR_FNR_calculations_MAF_threshold_local_concordance.r
```

### Input
This script requires `.validation` files, which are generated by:
```
validation_genotypes/validate_imputed_downsampled_vcfs.r
```

---

## Step 2: Calculate Heterozygous FPR in Sliding Windows

Sliding window analysis of heterozygous FPR is performed using the following script:

### Script
```
Local_heterozygous_FPR_calculations.sh
```

This script divides the genome into windows and calculates the heterozygous FPR for each.

---

## Step 3: Calculate Concordance in Sliding Windows

The concordance values across sliding windows are computed using:

### Script
```
Local_concordance_calculations.sh
```

This script:
- Calculates non reference concordance per window
- Identifies regions with low concordance and high FPR, which should be excluded from downstream analyses

---

## Step 4: Plotting Sliding Window Concordance and FPR

After generating individual output files for each sample, these can be visualized using Manhattan-style plots.

This is done with the individual output files from Local_heterozygous_FPR_calculations.sh and Local_concordance_calculations.sh
Files are prefixed with `Merged` (e.g., `Merged_acem2_x_x_x.bed`).

### Plotting Scripts Directory
```
plotting_scripts/
```

This directory contains R scripts for generating Manhattan plots of sliding window FPR and concordance.

---

## Notes
- Ensure all `.validation` files are available before beginning sliding window calculations.
- Scripts assume a consistent format for filenames and headers; edit as needed for custom data.
- Plotting scripts expect individual BED files with Merged prefix.


