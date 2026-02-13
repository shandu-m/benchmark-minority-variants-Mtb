#!/usr/bin/bash

if [ ${#@} -lt 6 ]; then
    echo "Usage: bash RUN_BATCH_get_bam_depths.sh [source_dir] [depth_output_dir] [scratch_dir] [hours] [memory] [config_file]"
    exit 1
fi

source_dir=$1
depth_output_dir=$2; mkdir -p $depth_output_dir

scratch_dir=$3
job_scratch_dir=$scratch_dir/job_files
err_out_dir=$scratch_dir/err_out
mkdir -p $job_scratch_dir $err_out_dir

hours=$4; memory=$5

config_file=$6
config_list=()
while read -r config
do
    IFS=": " read -ra entries <<< $config
    config_list+=(${entries[1]})

done < $config_file
script_dir=${config_list[0]}
email=${config_list[1]}
mail_type=${config_list[2]}

for isolate_tag_map in $(ls $source_dir)
do

    job_prefix=bam_depth_${isolate_tag_map}

    job_file=$job_scratch_dir/$job_prefix.sh
    > $job_file

    echo "#!/bin/bash" >> $job_file
    echo "" >> $job_file
    echo "#SBATCH -c 1" >> $job_file
    echo "#SBATCH -t 0-$hours:00" >> $job_file
    echo "#SBATCH -p short" >> $job_file
    echo "#SBATCH --mem=${memory}" >> $job_file
    echo "#SBATCH --mail-type=$mail_type" >> $job_file
    echo "#SBATCH --mail-user=$email" >> $job_file
    echo "#SBATCH -o $err_out_dir/${job_prefix}_%j.out" >> $job_file
    echo "#SBATCH -e $err_out_dir/${job_prefix}_%j.err" >> $job_file
    echo "" >> $job_file
    echo "module load gcc/14.2.0 samtools/1.21" >> $job_file
    echo "" >> $job_file

    while read -r tag_entry
    do
        IFS="," read -ra tag_entry_fields <<< $tag_entry
        sample_ID=${tag_entry_fields[0]}
        bam_filepath=${tag_entry_fields[1]} #merged BAM

        depth_file=$depth_output_dir/$sample_ID.depth.tsv
        depth_file_gzip=$depth_output_dir/$sample_ID.depth.tsv.gz

        echo "#$sample_ID" >> $job_file

        echo "samtools depth -a -Q 1 $bam_filepath > $depth_file" >> $job_file

        echo "gzip -c $depth_file > $depth_file_gzip" >> $job_file

        echo "rm $depth_file" >> $job_file
        
        echo "" >> $job_file

    done < $source_dir/$isolate_tag_map

    sbatch $job_file

done
