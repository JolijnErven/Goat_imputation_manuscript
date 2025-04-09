
awk -f get_size_bcftools.awk $1 > size_$1

sed -i 's/8.0-16.0_Mb/#C21807/g' size_$1
sed -i 's/4.0-8.0_Mb/#FF817E/g' size_$1
sed -i 's/2.0-4.0_Mb/#D1EAF0/g' size_$1
sed -i 's/1.0-2.0_Mb/#0077B6/g' size_$1
sed -i 's/0.5-1.0_Mb/#012A4A/g' size_$1
sed -i 's/>16_Mb/#800000/g' size_$1

cat head_bcftools.txt size_$1 > head_size_$1

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' head_size_$1 > final_size_$1

rm size_$1
rm head_size_$1
