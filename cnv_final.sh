#!/bin/bash

for i in *repeat_sanger.xls; 

do

echo "INFO: Processing $i"

name=${i%_ME*}

awk -v name=$name '{OFS="\t"}BEGIN{print name}NR>1{print $3}' $i > ${i%_ME*}_max_cov.tsv ;

done

region=$(ls -t *repeat_sanger.xls | head -1)

awk 'BEGIN{print "region"}NR>1{print $1}' $region | tr "," "\t" | awk '{print $1}' > region_tmp.tsv

paste region_tmp.tsv *max_cov.tsv > all_samples_max_cov.tsv


awk -v name=$name 'function divn(n, i, fn) {
    fn=n ".temp"
    printf "%s", $1 > fn
    for(i=2; i<=NF; i++)
       if (i != n)
          printf "%s%s", OFS, ($i ~ /[^0-9]/ ? $n "/" $i : ($i==0 ? "NA" : sprintf("%.3f", $n/$i) )) > fn
    print "" > fn
}
BEGIN {
   FS=OFS="\t"
}
{
   for (i=2; i<=NF; i++)
      divn(i)
}' all_samples_max_cov.tsv


for i in *.temp;

do 

echo "INFO: Renaming: $i"

name=$(head -1 $i |  awk '{for(i=3;i<=NF;++i)print $i}' | awk '{split($1,arr1,"/"); print arr1[1]}' | sort | uniq)

mv "$i" "$name"_final.tsv

done

