🐐 Ancient Goat Genomes – Imputation and Downstream Analyses
This repository contains all scripts and workflows used in the analysis of imputed ancient goat genomes as presented in the associated research paper.

📘 Erven, Jolijn A. M., Alice Etourneau, Marjan Mashkour, Mahesh Neupane, Phillipe Bardou, Alessandra Stella, Andrea Talenti, Clet Wandui Masiga, Curt Van Tassell, Emily Clark, François Pompanon, Licia Colli, Marcel Amillis, Marco Milanesi, Paola Crepaldi, The VarGoats Consortium, Bertrand Servin, Ben Rosen, Gwenola Tosser-Klopp, and Kevin G. Daly. 2025. “Inferring Domestic Goat Demographic History through Ancient Genome Imputation.” bioRxiv. doi:10.1101/2025.04.18.649576.

📂 Repository Structure & Workflow
Each subdirectory contains specific processing or analysis steps, from raw VCF handling to ROH, IBD, and FPR/concordance analyses.
There are prerequisites for some subdirectories, but these are mentioned in their README.md

## 🧑‍🔬 Workflow Overview
1. **`downsampling_BAM/`**: Downsample BAM files to various coverages.
2. **`Variant_calling/`**: Perform variant calling for high and low-coverage samples.
3. **`GLIMPSE2/`**: Genotype imputation using **GLIMPSE2**.
4. **`Validation_genotypes/`**: Validate imputed genotypes by concordance comparisons.
5. **`Genome-wide_concordance_FPR/`**: Sliding window analyses for concordance and FPR.
6. **`PCA/`**: Principal Component Analysis (PCA) of imputed genotypes.
7. **`F3/`**: F3 outgroup statistics to infer population relationships.
8. **`ROH/`**: Calculate **ROH** (Runs of Homozygosity) with **Plink** and **bcftools**.
9. **`local_ROH_heterozygosity/`**: Create local ROH and heterozygosity profiles.
10. **`IBD/`**: Perform **Refined IBD** analysis to classify individual relationships.

🧪 Reproducibility
All scripts are intended to be run on a Linux-based HPC environment using standard bioinformatics tools (GLIMPSE2, bcftools, plink, bcftools roh, bedtools, R, etc). Parameters and paths will need to be adjusted per user's dataset and system.

📄 Citation
If you use this repository or self-made scripts in your own research, please cite the corresponding publication:
📘 Erven, Jolijn A. M., Alice Etourneau, Marjan Mashkour, Mahesh Neupane, Phillipe Bardou, Alessandra Stella, Andrea Talenti, Clet Wandui Masiga, Curt Van Tassell, Emily Clark, François Pompanon, Licia Colli, Marcel Amillis, Marco Milanesi, Paola Crepaldi, The VarGoats Consortium, Bertrand Servin, Ben Rosen, Gwenola Tosser-Klopp, and Kevin G. Daly. 2025. “Inferring Domestic Goat Demographic History through Ancient Genome Imputation.” bioRxiv. doi:10.1101/2025.04.18.649576.

