#!/bin/bash

if [ ${#@} -lt 3 ]
then
    printf "\033[1mUsage: bash varscan.sh [tag_map] [work_dir] [pvalue]\033[0m"
    exit 1
fi

tag_map=$1
work_dir=$2
pvalue=$3

while read -r tag_entry
do
    IFS=',' read -ra tag_entry_fields <<< $tag_entry
    tag=${tag_entry_fields[0]}
    pileup_filepath=${tag_entry_fields[1]}

    printf "\n$tag: running VarScan...\n"

    out_snp_file=$work_dir/${tag}.snp.varscan.tsv
    out_indel_file=$work_dir/${tag}.indel.varscan.tsv

    java -jar /home/sm624/bin/VarScan.v2.3.9.jar mpileup2snp \
        $pileup_filepath \
        --min-var-freq 0.01 \
        --p-value $pvalue \
        > $out_snp_file

    java -jar /home/sm624/bin/VarScan.v2.3.9.jar mpileup2indel \
        $pileup_filepath \
        --min-var-freq 0.01 \
        --p-value $pvalue \
        > $out_indel_file

    out_file=$work_dir/${tag}.varscan.tsv
    cat $out_snp_file > $out_file
    awk 'FNR>1' $out_indel_file >> $out_file
    rm $out_snp_file $out_indel_file

done < $tag_map

printf "\nDONE\n\n"
