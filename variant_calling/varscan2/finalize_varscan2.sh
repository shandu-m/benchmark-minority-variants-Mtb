#!/usr/bin/bash

if [ ${#@} -lt 5 ]; then
    echo "Usage: bash finalize_varscan2.sh [source_dir] [output_dir] [scratch_dir] [hours] [memory]"
    exit 1
fi

source_dir=$1
output_dir=$2; mkdir -p $output_dir

scratch_dir=$3
job_scratch_dir=$scratch_dir/job_files
err_out_dir=$scratch_dir/err_out
temp_dir=$scratch_dir/temp
mkdir -p $job_scratch_dir $err_out_dir $temp_dir

hours=$4
memory=$5

tool="varscan2"

for isolate_tag_map in $(ls $source_dir)
do

    job_prefix=VS_finalize_${isolate_tag_map}

    job_file=$job_scratch_dir/$job_prefix.sh
    > $job_file

    echo "#!/bin/bash" >> $job_file
    echo "" >> $job_file
    echo "#SBATCH -c 1" >> $job_file
    echo "#SBATCH -t 0-$hours:00" >> $job_file
    echo "#SBATCH -p short" >> $job_file
    echo "#SBATCH --mem=${memory}" >> $job_file
    echo "#SBATCH --mail-type=FAIL" >> $job_file
    echo "#SBATCH --mail-user=smulaudzi@g.harvard.edu" >> $job_file
    echo "#SBATCH -o $err_out_dir/${job_prefix}_%j.out" >> $job_file
    echo "#SBATCH -e $err_out_dir/${job_prefix}_%j.err" >> $job_file
    echo "" >> $job_file

    while read -r tag_entry
    do
        IFS="," read -ra tag_entry_fields <<< $tag_entry
        tag=${tag_entry_fields[0]}
        raw_variants_file=${tag_entry_fields[1]}

        final_variants_file=${output_dir}/${tag}.${tool}.csv

        echo "#$sample_ID" >> $job_file

        ### STEP 1: finalize variants list ###

        echo "python /home/sm624/projects/mixed_calls/benchmarking/varscan/scripts/finalize_varscan2.py -v ${raw_variants_file} -o ${final_variants_file}" >> $job_file

        echo "" >> $job_file

    done < $source_dir/$isolate_tag_map

    sbatch $job_file

done