# Benchmarking within-sample minority variant detection with short-read sequencing in *M. tuberculosis*

Mulaudzi, S., Kulkarni, S., Marin, M. G. & Farhat, M. R. Benchmarking within-sample minority variant detection with short-read sequencing in *M. tuberculosis*. Preprint at [https://doi.org/10.64898/2026.02.13.704885](https://doi.org/10.64898/2026.02.13.704885) (2026).

**This is the project repo describing all sequence data simulation, sequence processing, variant calling and data analysis associated with our publication.**

- `data_simulation`: Describes our data simulation process.
- `wgs_processing`: Describes how we processed the simulated sequencing data.
- `variant_calling`: Describes how the benchmarked variant callers were run on the simulated strains.
- `variant_filtering`: Describes how we used an error model and allele fraction adjustment to filter FreeBayes SNVs and INDEL variants respectively. The error model is provided here to be used on FreeBayes variant calls.
- `analysis`: All benchmarking analysis. The comprehensive low mappability regions for variant filtering can be found at `analysis/data/L1234.H37Rv.pupmapper.indelHomologRegions.RLC.geneHomologs.first500bp.bed`.
