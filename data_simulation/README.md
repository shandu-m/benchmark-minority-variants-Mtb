# Data simulation

Reference genomes in which variants were introduced are located in `./reference_genomes`.

## H37Rv modification
Modification of the H37Rv reference genome was performed in the `h37rv_modification.ipynb` notebook.

For each mutant haplotype (n=10), the H37Rv reference genome was modified to create a mutant reference fasta for data simulation. Summaries of the variants introduced in each of these mutants are located in `./mutant_summaries/H37Rv`.

Candidate variant positions were chosen from the following files:
- Drug resistance variants: `./data/WHO_drug_SNP_mapping.csv`
- Homopolymer tract variants: `./data/DR_homopolymers.csv`
- Low mappability variants: `./data/low_mappability_regions.csv`

## L1-4 modification
Modification of a non-H37Rv reference genome from lineages 1-4 was performed in the `L1234_modification.ipynb` notebook.

For each mutant, an L1/L2/L3/L4 reference genome was modified to create a mutant reference fasta for data simulation. A summary of the variants introduced in each of these mutants is located in `./mutant_summaries/L1234`.

Candidate drug resistance variant positions were chosen from the `./data/WHO_drug_SNP_mapping.csv` file.

## InSilicoSeq (ISS)
Sequence data was simulated with InSilicoSeq using the script `./scripts/ISS_simulate_data.sh`.

InSilicoSeq can simulate sequencing data from multiple input reference fastas, with different specified abundances for these input references.

Main command for simulating the data:
```bash
iss generate --cpus 20 \
  --genomes $multifasta_file \
  --model $error_model \
  --output $output_file \
  --n_reads $n_reads \
  --abundance_file $abundance_dir/abundance_file_${freq}_${mutant}
```

Parameters used:
- `error_model` = hiseq
- `n_reads` (number of reads):
  - 50x: 2091299
  - 100x: 4182599
  - 200x: 8365198
  - 400x: 16730395
  - 700x: 29278191
 
The abundance files contain two lines, each with a SeqID followed by a frequency to specify the abundance of this reference fasta in the simulated sequencing reads. For example, the abundance file for the H37Rv mutant 1 strain with variant frequency 10% would look like this:
```bash
NC_000962.3 0.9
NC_000962.3-mutant-1 0.1
```

The associated multi-fasta file for this example strain would contain the fasta sequences for NC_000962.3 (non-mutant reference fasta) and NC_000962.3-mutant-1 (mutant reference fasta).

Initially we simulated only a subset of the sequencing depth and variant frequency combinations, before adding additional combinations to complete *all* sequencing depth (50-700x) and variant frequency (1-50%) pairings with `./scripts/ISS_simulate_data_additional.sh`.

## ART
Sequence data was simulated with ART using the script `./scripts/ART_simulate_data.sh`.

ART simulates from one reference fasta. To create intermediate variant allele frequencies, we simulated X reads from a non-mutant reference fasta and Y reads from mutant reference fasta, combining these reads to reach the appropriate sequencing depth of this simulated strain.

Main commands for simulating the data:
```bash
# simulate non-mutant reads
art_illumina -p \
  -i $ref_fasta \
  -l 150 \
  -f $novar_n_reads \
  -m 200 -s 10 \
  -o $novar_output_file

# simulate mutant reads
art_illumina -p \
  -i $fasta_file \
  -l 150 \
  -f $var_n_reads \
  -m 200 -s 10 \
  -o $var_output_file

# combine the simulated fastq files
cat ${novar_output_file}1.fq ${var_output_file}1.fq > $R1_file
cat ${novar_output_file}2.fq ${var_output_file}2.fq > $R2_file
```

Parameters used:
- `novar_n_reads` and `var_n_reads` (fold coverage):
  - These were calculated with the `./scripts/get_read_num_parameters.py` script to ensure that the sum of non-mutant simulated reads and mutant simulated reads amounted to the number of reads required for the sequencing depth of the simulated strain.
  - For example: for a simulated strain with depth 200x and variant frequency 10%, `novar_n_reads=211` and `var_n_reads=23`.
  - These fold coverages were calculated in the `./frequency_combinations/generate_read_num_parameters.ipynb` notebook.
 
The `ref_fasta` input fasta file is the non-mutant fasta file and the `fasta_file` input fasta file is the mutant fasta file.

Initially we simulated only a subset of the sequencing depth and variant frequency combinations, before adding additional combinations to complete *all* sequencing depth (50-700x) and variant frequency (1-50%) pairings with `./scripts/ART_simulate_data_additional.sh`.
