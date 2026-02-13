# Variant calling

Variants in all H37Rv and L1-4 strains were called by each benchmarked variant caller in the same way. The clinical isolates used for the "Real variant analysis" were also processed in the same way. The Mtb reference genome (H37Rv) was used for all tools.

Below is a summary of the commands used to run each variant caller. The inputs are:
- `ref_fasta`: H37Rv reference genome
- `bam_filepath`: BAM file for a strain created with our WGS processing pipeline
- `pileup_filepath`: mpileup file for a strain created with our WGS processing pipeline

## BinoSNP
```bash
# run BinoSNP
binoSNP \
  -r $ref_fasta \
  -i $interval_file \
  -m $mutation_file \
  -o $out_dir \
  $bam_filepath
```
The interval and mutation files (`$interval_file` and `$mutation_file`) required by BinoSNP were created according to [these guidelines](https://github.com/ngs-fzb/binoSNP).

## FreeBayes
```bash
# run FreeBayes
freebayes \
  -f $ref_fasta \
  -p 1 --min-alternate-count 2 --min-alternate-fraction 0.01 \
  --min-mapping-quality 30 --min-base-quality 30 \
  -b $bam_filepath \
  -v $out_file
```

## LoFreq
```bash
# insert INDEL qualities into BAM (in order to call INDELs)
lofreq indelqual \
  --dindel \
  --ref $ref_fasta \
  --out $indel_qual_bam \
  $bam_filepath

# index BAM with inserted INDEL qualities
samtools index $indel_qual_bam

# run lofreq
lofreq call-parallel \
  --pp-threads 20 \
  --call-indels \
  --force-overwrite \
  --ref $ref_fasta \
  --out $out_file \
  $indel_qual_bam
```

## Mutect2
```bash
# run Mutect2
gatk --java-options "-Dsamjdk.compression_level=5 -Xms2G -Xmx4G" \
  Mutect2 \
  --force-active \
  -R $ref_fasta \
  -I $bam_filepath \
  --f1r2-tar-gz $tag_dir/${tag}.f1r2.tar.gz \
  -O $tag_dir/${tag}.unfiltered.vcf \
  --annotation StrandBiasBySample \
  --max-reads-per-alignment-start 75

# learn read orientation model
gatk --java-options "-Dsamjdk.compression_level=5 -Xms2G -Xmx4G" \
  LearnReadOrientationModel \
  -I $tag_dir/${tag}.f1r2.tar.gz \
  -O $tag_dir/${tag}.read_orientation_model.tar.gz

# filter Mutect2 calls for microbes
gatk --java-options "-Dsamjdk.compression_level=5 -Xms2G -Xmx4G" \
  FilterMutectCalls \
  -R $ref_fasta \
  -V $tag_dir/${tag}.unfiltered.vcf \
  --stats $tag_dir/${tag}.unfiltered.vcf.stats \
  --orientation-bias-artifact-priors $tag_dir/${tag}.read_orientation_model.tar.gz \
  -O $tag_dir/${tag}.mutect2.microbial_mode.vcf \
  --microbial-mode true
```

## Pilon
```bash
# run Pilon
java -Xmx16G -jar /home/sm624/bin/pilon-1.24.jar \
  --genome $ref_fasta \
  --bam $bam_filepath \
  --output $output_file_prefix \
  --variant \
  --minmq 40 \
  --mindepth 5
```

## VarDict
```bash
# run VarDict
VarDict -th 12 \
  -G $ref_fasta \
  -b $bam_filepath \
  -R Chromosome:1-4411532 | teststrandbias.R | var2vcf_valid.pl -S > $out_file
```

## VarScan2
```bash
# run VarScan2 to call SNPs
java -jar /home/sm624/bin/VarScan.v2.3.9.jar mpileup2snp \
  $pileup_filepath \
  --min-var-freq 0.01 \
  --p-value 0.01 \
  > $out_snp_file

# run VarScan2 to call INDELs
java -jar /home/sm624/bin/VarScan.v2.3.9.jar mpileup2indel \
  $pileup_filepath \
  --min-var-freq 0.01 \
  --p-value 0.01 \
  > $out_indel_file

# combine SNP and INDEL variant files
cat $out_snp_file > $out_file
awk 'FNR>1' $out_indel_file >> $out_file
```

## Variant set normalization and filtering
Steps taken to standardize variant output across tools can be found in the file with prefix filter/finalize in each tool's directory.

Variant output across tools was combined into a single CSV file for each simulated strain in the `make_variant_summary.py` script, and for each clinical isolate in the `make_real_variant_summary.py` script.
