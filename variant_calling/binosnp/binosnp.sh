#!/bin/bash

if [ ${#@} -lt 2 ]
then
    printf "\033[1mUsage: bash binosnp.sh [tag_map] [work_dir]\033[0m"
    exit 1
fi

tag_map=$1
work_dir=$2

while read -r tag_entry
do
    IFS=',' read -ra tag_entry_fields <<< $tag_entry
    tag=${tag_entry_fields[0]}
    bam_filepath=${tag_entry_fields[1]}
    interval_file=${tag_entry_fields[2]}
    mutation_file=${tag_entry_fields[3]}

    printf "\n$tag: running BinoSNP...\n"

    out_dir=$work_dir

    rm -rf $work_dir/${tag}* #otherwise binosnp will not overwrite old files

    /home/sm624/bin/binoSNP/binoSNP \
        -r /n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna \
        -i $interval_file \
        -m $mutation_file \
        -o $out_dir \
        $bam_filepath

done < $tag_map

printf "\nDONE\n\n"
