# Whole-genome sequence data processing

The fastq files for each simulated strain were processed as follows:

1. Read trimming from the 3' end to a quality of 25 (PRINSEQ-lite): `./scripts/trim_reads.sh`.
2. Alignment to the Mtb refrence genome H37Rv (BWA-MEM) and duplicate read removal with Picard: `./scripts/generate_bam.sh`.
3. mpileup file creation with SAMtools: `./scripts/mpileup.sh`.

The resulting BAM and mpileup files were used for variant calling.
