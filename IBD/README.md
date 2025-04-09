# Refined IBD Analysis Pipeline for Imputed Ancient Goat Genomes

This pipeline performs identity-by-descent (IBD) analysis on imputed ancient goat genomes using the Refined IBD algorithm. 

## Overview of Steps
1. **Merge Imputed VCFs**  
   Combines VCFs from published and 0.5X imputed ancient genomes into a single merged VCF file.
2. **Run Refined IBD (3 repetitions)**  
   Executes Refined IBD three times to mitigate phasing errors, using a custom wrapper script.
3. **Filter IBD Segments**  
   Retains IBD segments greater than 3 cM (later repeated for 4 cM).
4. **Extract Unique Pairs**  
   Identifies all unique individual pairs sharing IBD segments.
5. **Extract IBD Per Pair**  
   Extracts all IBD segments for each unique pair.
6. **Merge IBD to BED Format**  
   Merges overlapping IBD segments into BED files.
7. **Add Genetic Map Coordinates**  
   Annotates IBD segments with cM coordinates using a genetic map.
8. **Classify Relationships**  
   Classifies each individual pair as related or unrelated based on number and total length of IBD segments.

All steps are repeated for both 0.6 cM and 4 cM merged datasets.

## Dependencies
- `bcftools`
- `Refined IBD` (Beagle)
- `mergeBed` (bedtools)
- `parallel`
- Custom scripts: `run_refined_IBD_3_times.sh`, `cm_add.sh` found in `scripts`

## Inputs
- Imputed VCFs from `/Imputed/published/glimpse_Plink/` and `/Imputed/glimpse_Plink/`
- A list of sites: `combined_above_3x_glimpse2_GP99_MAF_0.05_TV_no_low_conc_no_single_chr18_2mb_removed.pos.txt`
- Genetic map files (for coordinate conversion)

## Outputs
- Filtered IBD segment files (>3 cM)
- Per-pair BED files of merged IBD segments
- Final relationship classification files with columns: individual1, individual2, segment count, total length (cM), classification (rel/un)

# The plotting is done on altered files which can be found in the Tables associated with the manuscript (Table S15 and S16)

