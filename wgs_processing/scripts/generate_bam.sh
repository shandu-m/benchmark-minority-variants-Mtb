#!/usr/bin/bash

if [ ${#@} -lt 4 ]
then
    printf "Usage: bash generate_bam.sh [tag] [fastq1] [fastq2] [output_dir]\n"
    exit 1
fi

tag=$1

fastq1=$2
fastq2=$3

output_dir=$4/bam
if [ ! -d $output_dir ]
then
    mkdir $output_dir
fi

output_prefix=$output_dir/$tag

module load gcc/6.2.0 bwa/0.7.17 samtools/1.15.1 picard/2.8.0

sam_file=${output_prefix}.sam
bwa mem -M -t 20 /n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna $fastq1 $fastq2 > $sam_file

sorted_bam=${output_prefix}.sorted.bam
java -Xmx16G -jar $PICARD/picard-2.8.0.jar SortSam INPUT=$sam_file OUTPUT=$sorted_bam SORT_ORDER=coordinate
rm $sam_file

sorted_duprem_bam=${output_prefix}.sorted.duprem.bam
java -jar $PICARD/picard-2.8.0.jar MarkDuplicates I=$sorted_bam O=$sorted_duprem_bam M=${output_prefix}.dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    ASSUME_SORT_ORDER=coordinate
rm $sorted_bam
samtools index $sorted_duprem_bam
