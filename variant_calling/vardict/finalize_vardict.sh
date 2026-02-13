#!/usr/bin/bash

if [ ${#@} -lt 5 ]; then
    echo "Usage: bash finalize_vardict.sh [source_dir] [output_dir] [scratch_dir] [hours] [memory] (chromosome:yes)"
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

chromosome=$6
if [ ${#chromosome} -gt 0 ]; then 
    chromosome_switch=true
else
    chromosome_switch=false
fi

tool="vardict"

for isolate_tag_map in $(ls $source_dir)
do

    job_prefix=VD_finalize_${isolate_tag_map}

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
    echo "module load gcc/14.2.0 bcftools/1.21" >> $job_file #for bcftools
    echo "" >> $job_file

    while read -r tag_entry
    do
        IFS="," read -ra tag_entry_fields <<< $tag_entry
        tag=${tag_entry_fields[0]}
        vcf_filepath=${tag_entry_fields[1]}

        echo "#$sample_ID" >> $job_file

        ### STEP 1: add contig information to VCF ###

        old_header_file=${temp_dir}/${tag}.old.hdr
        new_header_file=${temp_dir}/${tag}.new.hdr
        contig_vcf=${temp_dir}/${tag}.${tool}.contig.vcf

        if $chromosome_switch; then
            contig_file=/home/sm624/inventory/h37rv_chrom_contig
        else
            contig_file=/home/sm624/inventory/h37rv_contig
        fi
        echo "bcftools view -h ${vcf_filepath} > ${old_header_file}" >> $job_file
        echo "head -3 ${old_header_file} > ${new_header_file}" >> $job_file
        echo "cat ${contig_file} >> ${new_header_file}" >> $job_file
        echo "tail -n +4 ${old_header_file} >> ${new_header_file}" >> $job_file

        echo "bcftools reheader -h ${new_header_file} -o ${contig_vcf} ${vcf_filepath}" >> $job_file

        ### STEP 2: split any multi-allelic sites and normalize indels;
        # simplify complex alleles with vcfwave ###

        MAsplit_norm_file=${temp_dir}/${tag}.${tool}.MAsplit.norm.vcf
        vcfwave_MAsplit_norm_file=${temp_dir}/${tag}.${tool}.vcfwave.MAsplit.norm.vcf

        if $chromosome_switch; then
            ref_file=/n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.chromosome.fna
        else
            ref_file=/n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna
        fi
        echo "bcftools norm --multiallelics -any -f ${ref_file} ${contig_vcf} -Ov -o ${MAsplit_norm_file}" >> $job_file
        echo "vcfwave --quiet ${MAsplit_norm_file} > ${vcfwave_MAsplit_norm_file}" >> $job_file

        ### STEP 3: extract necessary info with SNPsift and write to tsv file ###

        initial_variants_file=${temp_dir}/${tag}.temp.${tool}.tsv
        final_variants_file=${output_dir}/${tag}.${tool}.tsv

        echo "SnpSift extractFields -e '' -s ';' ${vcfwave_MAsplit_norm_file} \
            POS REF ID ALT FILTER AF DP VD REFBIAS VARBIAS MQ ORIGIN > ${initial_variants_file}" >> $job_file
        
        echo "echo -e 'POS\tREF\tID\tALT\tFILTER\tAF\tDP\tAC\tREFBIAS\tVARBIAS\tMQ\tORIGIN' > ${final_variants_file}" >> $job_file
        echo "tail -n+2 ${initial_variants_file} >> ${final_variants_file}" >> $job_file

        ### STEP 4: finalize variants list ###
        
        echo "python /home/sm624/projects/mixed_calls/benchmarking/${tool}/scripts/finalize_${tool}.py -v ${final_variants_file}" >> $job_file

        echo "rm ${temp_dir}/${tag}* ${final_variants_file}" >> $job_file

        echo "" >> $job_file

    done < $source_dir/$isolate_tag_map

    sbatch $job_file

done