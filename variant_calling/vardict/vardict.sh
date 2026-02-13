#!/bin/bash

if [ ${#@} -lt 2 ]
then
    printf "\033[1mUsage: bash vardict.sh [tag_map] [work_dir]\033[0m"
    exit 1
fi

tag_map=$1
work_dir=$2

module load gcc/9.2.0 R/4.3.1

while read -r tag_entry
do
    IFS=',' read -ra tag_entry_fields <<< $tag_entry
    tag=${tag_entry_fields[0]}
    bam_filepath=${tag_entry_fields[1]}

    printf "\n$tag: running VarDict...\n"

    out_file=$work_dir/${tag}.vardict.vcf

    VarDict -th 12 \
        -G /n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna \
        -b $bam_filepath \
        -R Chromosome:1-4411532 | teststrandbias.R | var2vcf_valid.pl -S > $out_file

done < $tag_map

printf "\nDONE\n\n"
