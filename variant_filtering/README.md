# Variant filtering

This directory contains all data, scripts and notebooks used to:
1. Generate the read mapping and quality metrics used to describe the characteristics of false positive variant calls made by all tools.
2. Predict true variant probability of FreeBayes SNVs for SNV filtering.
3. Adjust the alelle fraction of FreeBayes INDELs for INDEL filtering.

## Read mapping and quality metrics

We obtained read mapping and quality metrics for all variant calls made by all tools.

Computed in the `./scripts/get_bam_depths.sh` and `./scripts/get_genome_metrics.sh` scripts:
- Coverage ratio: Ratio of total coverage at site to maximum of rolling average from left/right.
- Soft-clipped bases ratio: Number of soft-clipped bases at each position, normalized by coverage.
- Discordantly-aligned reads ratio: Number of discordantly paired reads (non-LR reads) at each site, normalized by coverage.

Computed in the `./scripts/get_pos_based_metrics.sh` script:
- Base quality: Average base quality of bases supporting the variant allele.
- Strand bias: deviation of the proportion of forward reads out of the total reads from 0.5 (for variant callers that compute this: FreeBayes, Mutect2, VarDict, VarScan2).

These read mapping and quality metrics were used to describe the characteristics of true versus false variants called by all tools and as input to the error model for FreeBayes SNV call filtering.

## Predict true variant probability of FreeBayes SNVs for SNV filtering with the error model

Steps are executed in `filter_SNP_variants.ipynb` and `benchmarking_filter_FB_SNP_variants.Rmd`. Model data is contained in `error_model/`.

For each strain:
1. Filter the set of FreeBayes variants to include only SNPs.
2. Extract the metrics used as fixed effects in the model for each model.
3. Fit the error model.

## Adjust the allele fraction of FreeBayes INDELs for INDEL filtering

Steps are executed in `filter_INDEL_variants.ipynb`.

For each strain:
1. Filter the set of FreeBayes variants to include only INDELs.
2. Compute INDEL-related metrics (`./scripts/get_indel_metrics.py`).
3. Adjust the allele fraction of these INDEL variants as described in the paper.

