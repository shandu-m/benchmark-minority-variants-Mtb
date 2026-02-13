#!/bin/bash

if [ ${#@} -lt 2 ]
then
    printf "\033[1mUsage: bash mutect2.sh [tag_map] [work_dir]\033[0m"
    exit 1
fi

tag_map=$1
work_dir=$2

module load gcc/9.2.0 samtools/1.15.1

RGLine="@RG\tID:Sample1\tSM:Sample1"

while read -r tag_entry
do
    IFS=',' read -ra tag_entry_fields <<< $tag_entry
    tag=${tag_entry_fields[0]}
    bam_filepath=${tag_entry_fields[1]}

    printf "\n$tag: running Mutect2...\n"

    ref_fasta=/n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna

    tag_dir=$work_dir/$tag; mkdir -p $tag_dir

    bam=$work_dir/${tag}.sorted.duprem.samplenames.bam

    samtools addreplacerg -r $RGLine -o $bam $bam_filepath
    samtools index $bam

    # Mutect2
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms2G -Xmx4G" \
        Mutect2 \
        --force-active \
        -R $ref_fasta \
        -I $bam \
        --f1r2-tar-gz $tag_dir/${tag}.f1r2.tar.gz \
        -O $tag_dir/${tag}.unfiltered.vcf \
        --annotation StrandBiasBySample \
        --max-reads-per-alignment-start 75

    # Learn read orientation model
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms2G -Xmx4G" \
        LearnReadOrientationModel \
        -I $tag_dir/${tag}.f1r2.tar.gz \
        -O $tag_dir/${tag}.read_orientation_model.tar.gz

    # Filter Mutect2 calls for microbes
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms2G -Xmx4G" \
        FilterMutectCalls \
        -R $ref_fasta \
        -V $tag_dir/${tag}.unfiltered.vcf \
        --stats $tag_dir/${tag}.unfiltered.vcf.stats \
        --orientation-bias-artifact-priors $tag_dir/${tag}.read_orientation_model.tar.gz \
        -O $tag_dir/${tag}.mutect2.microbial_mode.vcf \
        --microbial-mode true

    rm $bam ${bam}.bai

done < $tag_map

printf "\nDONE\n\n"
