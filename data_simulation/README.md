# Data simulation

## ISS
Sequence data was simulated with InSilicoSeq using the script `./scripts/ISS_simulate_data.sh`.

Main command for simulating the data:
```bash
iss generate --cpus 20 --genomes $multifasta_file \
  --model $error_model \
  --output $output_file \
  --n_reads $n_reads \
  --abundance_file $abundance_dir/abundance_file_${freq}_${mutant}
```
Parameters used:
- error_model = hiseq
- n_reads (number of reads):
  - 50x: 2091299
  - 100x: 4182599
  - 200x: 8365198
  - 400x: 16730395
  - 700x: 29278191

Initially we simulated only a subset of the sequencing depth and variant frequency combinations, before adding additional combinations to complete *all* sequencing depth (50-700x) and variant frequency (1-50%) pairings with `./scripts/ISS_simulate_data_additional.sh`.

## H37Rv simulation
