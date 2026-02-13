#!/usr/bin/bash

if [ ${#@} -lt 11 ]
then
    printf "Usage: bash simulate_data_additional.sh [multifasta_dir] [abundance_dir] [output_dir] [scratch_dir] [n_reads] [memory] [hours] [queue] [error_model] [coverage] [label]\n"
    printf "    error_model: 'hiseq' for regular model, or specify the path to a pre-built model\n"
    printf "    coverage: expected coverage e.g. 200 for 200x\n"
    exit 1
fi

IFS=" " read -ra freqs <<< $(cat /home/sm624/projects/mixed_calls/benchmarking/comprehensive/freqs)

multifasta_dir=$1
abundance_dir=$2
output_dir=$3; mkdir -p $output_dir
scratch_dir=$4; mkdir -p $scratch_dir
n_reads=$5
memory=$6
hours=$7
queue=$8
error_model=$9
coverage=${10}
label=${11}

for mutant in $(ls $multifasta_dir)
do
    multifasta_file=$multifasta_dir/$mutant/${mutant}_full.fasta

    for freq in ${freqs[@]}
    do

        wc_grep=$(grep "$coverage,$freq" /home/sm624/projects/mixed_calls/benchmarking/comprehensive/new_combinations.csv | wc -l | awk '{print $1}')
        if [ $wc_grep -eq 0 ]
        then
            continue
        fi

        job_file=$scratch_dir/iss_${coverage}_${freq}_${mutant}_${label}.sh
        > $job_file

        # add job file specific specifications
        echo "#!/bin/bash" >> $job_file
        echo "" >> $job_file
        echo "#SBATCH -c 20" >> $job_file
        echo "#SBATCH -t 0-${hours}:00" >> $job_file
        echo "#SBATCH -p $queue" >> $job_file
        echo "#SBATCH --mem=${memory}G" >> $job_file
        echo "#SBATCH --mail-type=FAIL" >> $job_file
        echo "#SBATCH --mail-user=smulaudzi@g.harvard.edu" >> $job_file
        echo "#SBATCH -o $scratch_dir/iss_${coverage}_${freq}_${mutant}_%j.out" >> $job_file
        echo "#SBATCH -e $scratch_dir/iss_${coverage}_${freq}_${mutant}_%j.err" >> $job_file
        echo "" >> $job_file

        mutant_dir=$output_dir/${mutant}_${coverage}_${freq}
        if [ ! -d $mutant_dir ]; then
            mkdir $mutant_dir; fi

        output_file=$mutant_dir/${mutant}_${coverage}_${freq} 

        echo echo Running InSilicoSeq... >> $job_file
        echo iss generate --cpus 20 --genomes $multifasta_file --model $error_model --output $output_file --n_reads $n_reads --abundance_file $abundance_dir/abundance_file_${freq}_${mutant} >> $job_file
        R1_file=${output_file}_R1.fastq
        R2_file=${output_file}_R2.fastq
        echo "pigz --force $R1_file $R2_file" >> $job_file
        echo "" >> $job_file

        sbatch $job_file

    done
done
