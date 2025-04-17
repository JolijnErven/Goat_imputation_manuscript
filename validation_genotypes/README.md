# VCF Validation and Concordance Pipeline

This directory contains scripts that process VCF files to calculate concordance, false positive/negative rates, and generate plots. The pipeline consists of three primary stages:

---

## Step 1: Validate Imputed and Downsampled VCFs

The first step is to transform VCF data into a format more suitable for analysis in R. This is performed using the `validate_imputed_downsampled_vcfs.r` script.

### Script
```
validate_imputed_downsampled_vcfs.r
```

### Usage
You can use the following shell loop to create commands for all chromosomes and samples:
The $TRUE parameter are the diploid calls created from the Variant_calling/run_variant_calling_filtering.sh

```bash
for CHR in {1..29}; do
  for TRUE in $(ls path_to_TRUTH_VCF/*chr${CHR}_*.vcf.gz); do
    SAMPLE=$(basename $TRUE | cut -f1 -d'_')
    for DOWNSAMPLED in $(ls ${SAMPLE}*chr${CHR}_*vcf | cut -f1 -d'.'); do
      SAMPLE=$(echo $DOWNSAMPLED | cut -f1 -d'_')
      COV=$(echo $DOWNSAMPLED | cut -f2 -d'_' | sed -e "s/-/./g" | sed -e "s/X//g")
      GP=$(echo $DOWNSAMPLED | cut -f3 -d'_')
      echo "Rscript validate_imputed_downsampled_vcfs.r $TRUE ${DOWNSAMPLED}.vcf ${SAMPLE} $COV $GP ${DOWNSAMPLED}.validation 2> ${DOWNSAMPLED}.error" >> parallel_validation.sh
    done
  done
done
```

Then run all jobs in parallel:

```bash
parallel -a parallel_validation.sh -j8
```

This generates `.validation` files used in the next stage.

---

## Step 2: Calculate Concordance, FPR, and FNR Metrics

The `.validation` files are processed using the following R scripts to compute concordance, false positive rate (FPR), false negative rate (FNR) for MAF tresholds and tranches

### Scripts
- `Concordance_FPR_FNR_calculations_MAF_threshold.r`
- `Concordance_FPR_FNR_calculations_MAF_tranches.r`

### Example Usage
```bash
Rscript Concordance_FPR_FNR_calculations_MAF_threshold.r output.validation
# and
Rscript Concordance_FPR_FNR_calculations_MAF_tranches.r output.validation
```

These scripts produce summary tables of validation metrics stratified by sample, coverage and GP for MAF thresholds and tranches.

---

## Step 3: Plotting

The output from Step 2 is visualized using the R scripts in the `plotting_scripts` directory.

There are two scripts in this directory
- Script to plot concordance, FPR and FNR (`Plot_validation_imputation.R`)
- Script to plot recovery (`Plot_imputation_recovery_MAF_threshold.R`)

## Notes
- Ensure all necessary R packages are installed before running the scripts.
- Scripts assume a specific filename structure for VCF files; you may need to adapt the parsing logic for different naming conventions.
- This pipeline is optimized for high-throughput validation across many samples and chromosomes.

