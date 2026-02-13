#!/usr/bin/bash

if [ ${#@} -lt 6 ]
then
    printf "Usage: bash trim_reads.sh [tag] [fastq1] [fastq2] [output_dir] [scratch_dir] [trim_type]\n"
    exit 1
fi

tag=$1

fastq1=$2
fastq2=$3

output_dir=$4/prinseq
if [ ! -d $output_dir ]
then
    mkdir $output_dir
fi

scratch_dir=$5

trim_type=$6

log_out=$scratch_dir/${tag}_prinseq.log

trimmed_out=$output_dir/${tag}.trimmed

if [ "$trim_type" == "mqm20" ]
then
    prinseq-lite.pl -fastq $fastq1 -fastq2 $fastq2 \
        -out_good $trimmed_out \
        -out_bad null \
        -log $log_out \
        --min_qual_mean 20 \
        -verbose
elif [ "$trim_type" == "mqm25" ]
then
    prinseq-lite.pl -fastq $fastq1 -fastq2 $fastq2 \
        -out_good $trimmed_out \
        -out_bad null \
        -log $log_out \
        --min_qual_mean 25 \
        -verbose
elif [ "$trim_type" == "mqm30" ]
then
    prinseq-lite.pl -fastq $fastq1 -fastq2 $fastq2 \
        -out_good $trimmed_out \
        -out_bad null \
        -log $log_out \
        --min_qual_mean 30 \
        -verbose
elif [ "$trim_type" == "tqr20" ]
then
    prinseq-lite.pl -fastq $fastq1 -fastq2 $fastq2 \
        -out_good $trimmed_out \
        -out_bad null \
        -log $log_out \
        --trim_qual_right 20 \
        -verbose
elif [ "$trim_type" == "tqr25" ]
then
    prinseq-lite.pl -fastq $fastq1 -fastq2 $fastq2 \
        -out_good $trimmed_out \
        -out_bad null \
        -log $log_out \
        --trim_qual_right 25 \
        -verbose
elif [ "$trim_type" == "tqr30" ]
then
    prinseq-lite.pl -fastq $fastq1 -fastq2 $fastq2 \
        -out_good $trimmed_out \
        -out_bad null \
        -log $log_out \
        --trim_qual_right 30 \
        -verbose
else
    echo "Trim option not supported!"
fi
