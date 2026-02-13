#!/bin/bash

if [ ${#@} -lt 2 ]
then
    printf "\033[1mUsage: bash freebayes.sh [tag_map] [work_dir]\033[0m"
    exit 1
fi

tag_map=$1
work_dir=$2

while read -r tag_entry
do
    IFS=',' read -ra tag_entry_fields <<< $tag_entry
    tag=${tag_entry_fields[0]}
    bam_filepath=${tag_entry_fields[1]}

    printf "\n$tag: running Freebayes...\n"

    out_file=$work_dir/${tag}.freebayes.vcf

    freebayes \
        -f /n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna \
        -p 1 --min-alternate-count 2 --min-alternate-fraction 0.01 \
        --min-mapping-quality 30 --min-base-quality 30 \
        -b $bam_filepath \
        -v $out_file

done < $tag_map

printf "\nDONE\n\n"
