#!/bin/bash


# PARAMETERS #  
thread=35  # number of threads for parallelisation 
CHROM=("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29") # list of chromosomes, autosomes only 
RUN=/Software/bcftools-1.17/    # Path to bcftools
dir=/imputation/  # working directory 
merged="Imputed"    # prefix
REF=//modern_reference_dataset/   #  path to reference panel in vcf.gz format
MAP=/vargoats/bertrand_servin_map_breed_sex/  #  path to genetic map

#pipeline followed from https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/
##############################################################################################################################

# STEP 1: SPLIT THE GENOME INTO CHUNKS
mkdir -p ${dir}/glimpse_ref_panel/
for CHR in ${CHROM}; do
if [ ! -f ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt ]; then
		/programs/GLIMPSE-2.0.0/GLIMPSE2_chunk_static  --sequential --input ${REF}chr${CHR}_$3.vcf.gz \
			--region ${CHR} --window-mb 4 --buffer-mb 0.2 --map ${MAP}_${CHR}.gmap \
			--output ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt --log  ${dir}/glimpse_ref_panel/chunks.chr${CHR}.log
	fi
done

# STEP 2: SPLIT THE GENOME INTO bins
for CHR in ${CHROM}; do

#	# Extract variable from chunk file
	while IFS="" read -r LINE || [ -n "$LINE" ];
         do
         printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
         IRG=$(echo $LINE | cut -d" " -f3)
         ORG=$(echo $LINE | cut -d" " -f4)

          /programs/GLIMPSE-2.0.0/GLIMPSE2_split_reference_static --reference ${REF}chr${CHR}_$3.vcf.gz --map ${MAP}_${CHR}.gmap \
		  --input-region ${IRG} --output-region ${ORG} --output ${dir}/glimpse_ref_panel/vargoats_snps_filtered_1372_beagle5-3_chr${CHR}_phased_glimpse2-filtering
     done < ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt
done

# STEP 3: IMPUTE AND PHASE
mkdir -p ${dir}/glimpse_Imputed

for CHR in ${CHROM}; do

#	# Extract variable from chunk file
	IDvar=`cut -f 1 ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt`
	REGS=`cut -f 3 ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt | cut -d":" -f 2 | cut -d"-" -f1`
	REGE=`cut -f 3 ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt | cut -d":" -f 2 | cut -d"-" -f2`

	# Convert string variable to list (for looping afterwards)
	IDlistA=($IDvar); IDlistB=($IDvar); REGSlist=($REGS); REGElist=($REGE);
	for ((i=0;i<10;i++)); do IDlistB[$i]="0"${IDlistB[$i]}; done # convert the single digits to double digits

	# Parallelise the corresponding number of threads
	# (If ID is not dividable by the number of threads or equal to zero, then run in background, else run in foregound:
	# this allows to run jobs on a controlled number threads.)
	for i in ${IDlistA[@]}; do
		if [ $(expr $i % $thread ) -ne 0 ] || [ $i -eq 0 ]; then
			/programs/GLIMPSE-2.0.0/GLIMPSE2_phase_static --bam-file ${dir}/chr${CHR}_$1.bam \
				--reference ${dir}/glimpse_ref_panel/vargoats_snps_filtered_1372_beagle5-3_chr${CHR}_phased_glimpse2-filtering_${REGSlist[$i]}_${REGElist[$i]}.bin \
				--output ${dir}/glimpse_Imputed/chr${CHR}_${REGSist[$i]}_${REGElist[$i]}.bcf &
  		# wait for the last chunk to finish before indexing
			if [ $i -eq ${IDlistA[-1]} ]; then
				wait
			fi
		else
          		/programs/GLIMPSE-2.0.0/GLIMPSE2_phase_static --bam-file ${dir}/chr${CHR}_$1.bam \
                     		--reference ${dir}/glimpse_ref_panel/vargoats_snps_filtered_1372_beagle5-3_chr${CHR}_phased_glimpse2-filtering_${REGSlist[$i]}_${REGElist[$i]}.bin \
                     		--output ${dir}/glimpse_Imputed/chr${CHR}_${REGSlist[$i]}_${REGElist[$i]}.bcf
		fi
	done
	
	# Index
	for i in ${IDlistA[@]}; do
		${RUN}bcftools index -f ${dir}/glimpse_Imputed/chr${CHR}_${REGSlist[$i]}_${REGElist[$i]}.bcf
	done
done

# STEP 4: LIGATE MULTIPLE CHUNKS
mkdir -p ${dir}/glimpse_Ligated
for CHR in ${CHROM}; do
	ls -1v ${dir}/glimpse_Imputed/chr${CHR}.*.bcf > ${dir}/glimpse_Ligated/list.chr${CHR}.txt
	/programs/GLIMPSE-2.0.0/GLIMPSE2_ligate_static --input ${dir}/glimpse_Ligated/list.chr${CHR}.txt \
		--output ${dir}/glimpse_Ligated/${merged}.chr${CHR}_$2.ligated.bcf
	${RUN}bcftools index -f ${dir}/glimpse_Ligated/${merged}.chr${CHR}_$2.ligated.bcf
done

# STEP 5: CONVERT TO VCF
for CHR in ${CHROM}; do
	# STEP 4.1: CONVERT BCF TO VCF
	${RUN}bcftools view -Ov \
		-o ${dir}/glimpse_Ligated/${merged}.chr${CHR}_$2.ligated.vcf \
		${dir}/glimpse_Ligated/${merged}.chr${CHR}_$2.ligated.bcf
	${RUN}bcftools index -f ${dir}/glimpse_Ligated/${merged}.chr${CHR}_$2.ligated.vcf

done

# STEP 6: REMOVE FILES
for CHR in ${CHROM}; do

	rm ${dir}/glimpse_Imputed/*
#	rm ${dir}/glimpse_ref_panel/*
#	rm ${dir}/glimpse_Ligated/${merged}.chr${CHR}.ligated.bcf*
done

echo "Your results consists of the following files:"
echo "${dir}/glimpse_Ligated/${merged}.chr${CHR}.ligated.vcf"
echo "END OF SCRIPT"



