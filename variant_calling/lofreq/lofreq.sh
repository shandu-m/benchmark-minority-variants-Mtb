#!/bin/bash

if [ ${#@} -lt 2 ]
then
    printf "\033[1mUsage: bash lofreq.sh [tag_map] [work_dir]\033[0m"
    exit 1
fi

tag_map=$1
work_dir=$2

ref_fasta=/n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna

while read -r tag_entry
do
    IFS=',' read -ra tag_entry_fields <<< $tag_entry
    tag=${tag_entry_fields[0]}
    bam_filepath=${tag_entry_fields[1]}

    printf "\n$tag: running LoFreq...\n"

    indel_qual_bam=$work_dir/${tag}.sorted.indelqual.bam; rm -f $indel_qual_bam
    out_file=$work_dir/${tag}.lofreq.vcf
    
    lofreq indelqual --dindel --ref $ref_fasta --out $indel_qual_bam $bam_filepath
    samtools index $indel_qual_bam
    lofreq call-parallel --pp-threads 20 --call-indels --force-overwrite --ref $ref_fasta --out $out_file $indel_qual_bam
    
    rm $indel_qual_bam ${indel_qual_bam}.bai

done < $tag_map

printf "\nDONE\n\n"
