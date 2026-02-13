#!/bin/bash

if [ ${#@} -lt 2 ]
then
    printf "\033[1mUsage: bash pilon.sh [tag_map] [work_dir]\033[0m"
    exit 1
fi

tag_map=$1
work_dir=$2

while read -r tag_entry
do
    IFS=',' read -ra tag_entry_fields <<< $tag_entry
    tag=${tag_entry_fields[0]}
    bam_filepath=${tag_entry_fields[1]}

    printf "\n$tag: running Pilon...\n"
    
    output_prefix=$work_dir/${tag}.pilon
    
    java -Xmx16G -jar /home/sm624/bin/pilon-1.24.jar \
        --genome /n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna \
        --bam $bam_filepath \
        --output $output_prefix \
        --variant \
        --minmq 40 \
        --mindepth 5 
    
    rm ${output_prefix}.fasta

done < $tag_map

printf "\nDONE\n\n"
