#!/usr/bin/bash

if [ ${#@} -lt 2 ]
then
    printf "\n\033[1mUsage: bash mpileup.sh [tag_path_map] [output_dir]\033[0m\n\n"
    exit 1
fi

tag_path_map=$1; output_dir=$2

module load gcc/9.2.0 samtools/1.15.1

reference=/n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna 

while read -r tag_entry
do
    IFS="," read -ra fields <<< $tag_entry
    tag=${fields[0]}
    bam=${fields[1]}
    
    echo ">>> $tag <<<"
    >&2 echo ">>> $tag <<<"
    
    pileup_file=$output_dir/${tag}.mpileup
    
    samtools mpileup -f $reference $bam -o $pileup_file
    
    echo ">>> $tag:DONE <<<"
    >&2 echo ">>> $tag:DONE <<<"

done < $tag_path_map
