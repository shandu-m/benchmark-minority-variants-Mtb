# Analysis

This directory contains all notebooks used for benchmarking analysis.

## Main results
- `0_simulated_data_statistics.ipynb`: analysis of the simulated data.
- `A_overall_accuracy.ipynb`: analysis of overall tool accuracy.
- `B_accuracy_by_AF_and_depth.ipynb`: analysis of tool accuracy by increasing variant frequency and depth.
- `C_FP_characteristics.ipynb`: analysis of false positive variant call characteristics.
- `D_FB_FP_filtering.ipynb`: analysis of different filtering schemes for the FreeBayes variants.
- `E_AF_accuracy.ipynb`: analysis of variant caller allele frequency accuracy.
- `F_binosnp_data_analysis.ipynb`: analysis of variant caller performance in *in-vitro* data.

Each of the above notebooks has a corresponding source notebook in which data tables etc. were created.

## Supplementary results
- `SUPP_A_binosnp_analysis.ipynb`: analysis of BinoSNP performance (Supplementary Results B).
- `SUPP_B_problematic_lineage_SNPs.ipynb`: analysis of variant caller failure to detect simulated variants in close proximity to fixed baseline lineage variants (Supplementary Results C).
- `SUPP_C_ART_comparison.ipynb`: analysis of variant caller performance in the ART-simulated data (Supplementary Results E).

Each of the above notebooks has a corrsponding source notebook in which data tables etc. were created.

## Methods
- `METHODS_AF0.ipynb`: analysis of variant call set difference between min AF = 0% and min AF = 1%.
- `METHODS_real_variant_analysis.ipynb`: analysis of variants called by all tools in a set of clinical isolates.

Analysis for Supplementary Results A is in `0_simulated_data_statistics.ipynb`. Analysis for Supplementary Results D is in `METHODS_real_variant_analysis.ipynb`.

The `data` folder contains files used for annotating genomic regions according to drug resistance (DR), homopolymer tract (HT) and low mappability (LM) regions. The `LM_pos_list.pkl` file corrresponds to the candidate positions used for LM variant simulation. The `L1234.H37Rv.pupmapper.indelHomologRegions.RLC.geneHomologs.first500bp.bed` file corresponds to the comprehensive low mappability regions used for false positive characterization and region masking for variant filtering.
