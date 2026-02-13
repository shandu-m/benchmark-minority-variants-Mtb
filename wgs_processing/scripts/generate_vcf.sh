#!/usr/bin/bash

if [ ${#@} -lt 3 ]
then
    printf "Usage: bash generate_vcf.sh [tag] [bam] [output_dir]\n"
    exit 1
fi

tag=$1

bam=$2

output_dir=$3/pilon
if [ ! -d $output_dir ]
then
    mkdir $output_dir
fi

output_prefix=$output_dir/$tag

java -Xmx16G -jar /home/sm624/bin/pilon-1.23.jar \
    --genome /n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna \
    --bam $bam \
    --output $output_prefix \
    --variant --threads 20

gzip -f ${output_prefix}.vcf
