# ROH Datasets and Analysis Pipeline

This directory contains all scripts used to generate Runs of Homozygosity (ROH) datasets for the study. 
Each script is stand-alone, meaning the order of execution does not matter. However, to produce the final plots, all relevant scripts should be executed beforehand.

## Overview
- All scripts are self-contained and can be run independently.
- Final plots require outputs from all analysis scripts.
- Scripts explore both PLINK and BCFtools-based ROH calling strategies under varying coverage and conditions.

---

## Script Summaries

### `Plink_ROH_pipeline_SNP_exploration.sh`
Explores the minimum number of SNPs required to call a ROH using PLINK. It also tests the effect of varying levels of downsampling.

- **Goal:** Assess a SNP density threshold impacts ROH detection accuracy.
- **Output:** Dataset for **Figure S28** in the study.

---

### `ROH_calcuation_individual_no_missingness.sh`
Calculates ROH using PLINK on individual test genomes imputed at different depths (0.25X, 0.5X, 0.75X, 1X) and the non imputed high coverage.

- **Goal:** Determine if imputation introduces false genotypes that alter ROH estimation.
- **Output:** Dataset for **Figure S29** in the study.

---

### `BCFtools_ROH_no_missingness_imputed_test_samples.sh`
Performs ROH calling with BCFtools on high coverage and 0.5X imputed test genomes.

- **Goal:** Evaluate imputation accuracy in ROH detection using BCFtools.
- **Output:** Dataset for **Figure S30** in the study.

---

### `ROH_ancient_imputed_paleogenomes_PLINK_BCFtools.sh`
Runs both PLINK and BCFtools to compute ROH on published paleogenomes (>0.5X coverage).

- **Goal:** Cross-compare ROH calls between tools on ancient samples to infer demographic history.
- **Output:** Dataset for **Figure S32** and a figure on OSF [10.17605/osf.io/av4f9](https://doi.org/10.17605/osf.io/av4f9).

---

### `Merge_4Mb_Bcftools_Plink_ROH.sh`
Combines large ROHs (>4 Mb) called by BCFtools into the PLINK ROH dataset for imputed paleogenomes.

- **Goal:** Merge high-confidence long ROHs from bcftools onto plink.
- **Output:** Dataset for **Figure 5** in the main paper.

---

## Plotting
Once all scripts have been executed, you can generate visualizations using the provided plotting script (in `plotting_scripts`). Make sure all intermediate ROH results are present before plotting.

---

## Notes
- Ensure all necessary tools (`plink`, `bcftools`, etc.) are installed.
- Input VCF and BIM/FAM/BED files should be formatted appropriately for each tool.


