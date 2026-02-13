#!/usr/bin/bash

if [ ${#@} -lt 8 ]
then
    printf "Usage: bash simulate_data.sh [fasta_dir] [output_dir] [scratch_dir] [memory] [hours] [queue] [coverage] [label] ..\n"
    printf "    fasta_dir: path to directory in which fasta files from which to simulate are stored\n"
    printf "    output_dir: path to directory in which to store the simulated read files\n"
    printf "    scratch_dir: path to directory in which to store temporary files\n"
    printf "    optional parameters:\n"
    printf "      --lineage: specify path to alternative lineage references\n"
    exit 1
fi

IFS=" " read -ra freqs <<< $(cat /home/sm624/projects/mixed_calls/benchmarking/comprehensive/freqs)

fasta_dir=$1
output_dir=$2; mkdir -p $output_dir
scratch_dir=$3; mkdir -p $scratch_dir
memory=$4
hours=$5
queue=$6
coverage=$7
label=$8

if [ ${#@} -gt 8 ] 
then
    lineage=$9; ref_dir=${10}
else
    ref_fasta=/n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna; fi
    
for mutant in $(ls $fasta_dir)
do
    fasta_file=$fasta_dir/$mutant/${mutant}.fasta

    if [ ${#@} -gt 8 ]; then
        IFS="_" read -ra mutant_fields <<< $mutant
        mutant_lineage=${mutant_fields[0]}
        ref_fasta=$ref_dir/$mutant/${mutant_lineage}.fasta
    fi

    for freq in ${freqs[@]}
    do

        novar_n_reads=$(python /home/sm624/projects/mixed_calls/benchmarking/comprehensive/get_read_num_parameters.py "ART" $coverage $freq "no_var")
        if [ "$novar_n_reads" == "" ]
        then
            continue
        else
            var_n_reads=$(python /home/sm624/projects/mixed_calls/benchmarking/comprehensive/get_read_num_parameters.py "ART" $coverage $freq "var")
        fi
        
        mutant_novar=$mutant"_novar"
        mutant_var=$mutant"_var"

        mutant_dir=$output_dir/${mutant}_${coverage}_${freq}
        if [ ! -d $mutant_dir ]; then
            mkdir $mutant_dir; fi

        novar_output_file=$mutant_dir/${mutant_novar}_${coverage}_${freq}_
        var_output_file=$mutant_dir/${mutant_var}_${coverage}_${freq}_

        R1_file=$mutant_dir/${mutant}_${coverage}_${freq}_R1.fastq
        R2_file=$mutant_dir/${mutant}_${coverage}_${freq}_R2.fastq

        job_file=$scratch_dir/art_${coverage}_${freq}_${mutant}_${label}.sh
        > $job_file

        echo "#!/bin/bash" >> $job_file
        echo "" >> $job_file
        echo "#SBATCH -c 1" >> $job_file
        echo "#SBATCH -t 0-${hours}:00" >> $job_file
        echo "#SBATCH -p $queue" >> $job_file
        echo "#SBATCH --mem=${memory}M" >> $job_file
        echo "#SBATCH --mail-type=FAIL" >> $job_file
        echo "#SBATCH --mail-user=smulaudzi@g.harvard.edu" >> $job_file
        echo "#SBATCH -o $scratch_dir/art_${coverage}_${freq}_${mutant}_%j.out" >> $job_file
        echo "#SBATCH -e $scratch_dir/art_${coverage}_${freq}_${mutant}_%j.err" >> $job_file
        echo "" >> $job_file

        echo "art_illumina -p -i $ref_fasta -l 150 -f $novar_n_reads -m 200 -s 10 -o $novar_output_file" >> $job_file

        echo "art_illumina -p -i $fasta_file -l 150 -f $var_n_reads -m 200 -s 10 -o $var_output_file" >> $job_file

        echo "cat ${novar_output_file}1.fq ${var_output_file}1.fq > $R1_file" >> $job_file
        echo "cat ${novar_output_file}2.fq ${var_output_file}2.fq > $R2_file" >> $job_file

        echo "rm ${novar_output_file}1.fq ${novar_output_file}1.aln ${novar_output_file}2.fq ${novar_output_file}2.aln" >> $job_file
        echo "rm ${var_output_file}1.fq ${var_output_file}1.aln ${var_output_file}2.fq ${var_output_file}2.aln" >> $job_file

        echo "pigz --force $R1_file $R2_file" >> $job_file

        sbatch $job_file

    done

done
