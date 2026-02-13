#!/home/sm624/.conda/envs/python/bin/python

import pandas as pd
import numpy as np
import os
import sys
import argparse, warnings

def finalize_lofreq(variants_file):

    variants_df = pd.read_csv(variants_file, sep="\t")

    ### filter for AC ≥ 2 and AF ≥ 0.01 ###
    
    DPs = variants_df.DP.values
    AFs = variants_df.AF.values

    ACs = [AF*DP for AF,DP in zip(AFs, DPs)]
    variants_df["AC"] = ACs

    filtered_variants_df = variants_df[(variants_df.AC >= 2) & (variants_df.AF >= 0.01)]

    filtered_variants_df = filtered_variants_df.sort_values("POS")
    
    variants_file_out = variants_file[:-3] + "csv"
    variants_df.to_csv(variants_file_out, index=False)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-v", "--variantFile", dest='variants_file', type=str, required=True, help='TSV file with variants pulled from a normalized LoFreq VCF to be finalized.')
    
    cmd_line_args = parser.parse_args()
    
    variants_file = cmd_line_args.variants_file

    finalize_lofreq(variants_file)