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


awk '
{
   gsub(/\r/,"")
}
{OFS="\t"}{
  nf=NF
  close(out_file)
  for(k=2;k<=nf;k++){
    out_file=""
    for(i=2;i<=nf;i++){
      if($i!=0){
         $(NF+1)=sprintf("%.03f",$k/$i)
      }
      else{
         $(NF+1)=sprintf("%s","NaN")
      }
    }
    out_file=k"field_out_file.tsv"
    print >> (out_file)
    NF=nf
  }
}' all_samples_max_cov.tsv
