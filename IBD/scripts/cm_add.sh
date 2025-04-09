#!/bin/bash

# This script takes an IBD file (5 columns: ind1 ind2 chr pos1 pos2) and:
# For each row, locates the closest matching position in a genetic map.
# Converts physical positions (pos1, pos2) to genetic positions (cM).
# Outputs the new data to input_file.cm.
# Script created by Lara Cassidy and modified by Kevin Daly

file=${1?Error: Give input *.ibd file - five columns (pair1,pair2,chr,pos1,pos2)}
for line in $(cat $file | awk '{print $1"+"$2"+"$3"+"$4"+"$5}'); do chr=$(echo $line | cut -f3 -d+); start=$(echo $line | cut -f4 -d+); end=$(echo $line | cut -f5 -d+); ind1=$(echo $line | cut -f1 -d+); ind2=$(echo $line | cut -f2 -d+)
mb_start=$(cat /bertrand_servin_map_breed_sex/Recombination_maps_1Mb_global_population_breed_sex-averaged_chr${chr}.gmap | awk -v t=$start -v c=4  '{a[NR]=$c}END{asort(a);d=a[NR]-t;d=d<0?-d:d;v = a[NR]; for(i=NR-1;i>=1;i--){ m=a[i]-t;m=m<0?-m:m; if(m<d){ d=m;v=a[i] } } print v }')
mb_end=$(cat bertrand_servin_map_breed_sex/Recombination_maps_1Mb_global_population_breed_sex-averaged_chr${chr}.gmap | awk -v t=$end -v c=4  '{a[NR]=$c}END{asort(a);d=a[NR]-t;d=d<0?-d:d;v = a[NR]; for(i=NR-1;i>=1;i--){ m=a[i]-t;m=m<0?-m:m; if(m<d){ d=m;v=a[i] } } print v }')
cm_start=$(cat bertrand_servin_map_breed_sex/Recombination_maps_1Mb_global_population_breed_sex-averaged_chr${chr}.gmap | awk -v var=$mb_start '{if ($4 == var) print $3}')
cm_end=$(cat bertrand_servin_map_breed_sex/Recombination_maps_1Mb_global_population_breed_sex-averaged_chr${chr}.gmap | awk -v var=$mb_end '{if ($4 == var) print $3}')
echo $ind1 $ind2 $chr $start $end $cm_start $cm_end | awk '{print $0}' ; done >$file.cm
