#!/home/sm624/.conda/envs/python/bin/python

import pandas as pd
import numpy as np
import os
import sys
import argparse, warnings

def finalize_varscan2(variants_file, output_file):

    raw_variants_df = pd.read_csv(variants_file, sep="\t")

    ### reformat ###
    variants_df = pd.DataFrame(columns=["POS", "REF", "ALT", "AF", "DP", "AO_F", "AO_R"])

    variants_df["POS"] = raw_variants_df.Position
    variants_df["REF"] = raw_variants_df.Ref
    variants_df["ALT"] = raw_variants_df.Var
    variants_df["AF"] = [float(cov_field.split(":")[-2].strip("%"))/100 for cov_field in raw_variants_df["Cons:Cov:Reads1:Reads2:Freq:P-value"]]
    variants_df["DP"] = [int(cov_field.split(":")[1]) for cov_field in raw_variants_df["Cons:Cov:Reads1:Reads2:Freq:P-value"]]
    variants_df["AO_F"] = [int(read_field.split(":")[-3]) for read_field in raw_variants_df["StrandFilter:R1+:R1-:R2+:R2-:pval"]]
    variants_df["AO_R"] = [int(read_field.split(":")[-2]) for read_field in raw_variants_df["StrandFilter:R1+:R1-:R2+:R2-:pval"]]
    
    ### filter for AC ≥ 2 ###
    DPs = variants_df.DP.values
    AFs = variants_df.AF.values

    ACs = [AF*DP for AF,DP in zip(AFs, DPs)]
    variants_df["AC"] = ACs

    variants_df = variants_df[variants_df.AC >= 2]

    ### handle INDELs ###

    ins_df = variants_df[variants_df.ALT.str.contains("+", regex=False)].copy(deep=True).reset_index(drop=True)
    ins_REFs = ins_df.REF.values
    ins_ALTs = ins_df.ALT.values
    new_ins_ALTs = [ref + alt.strip("+") for ref, alt in zip(ins_REFs, ins_ALTs)]
    ins_df["ALT"] = new_ins_ALTs
    
    del_df = variants_df[variants_df.ALT.str.contains("-", regex=False)].copy(deep=True).reset_index(drop=True)
    del_REFs = del_df.REF.values
    del_ALTs = del_df.ALT.values
    new_del_REFs = [ref + alt.strip("-") for ref, alt in zip(del_REFs, del_ALTs)]
    new_del_ALTs = [ref[0] for ref in new_del_REFs]
    del_df["REF"] = new_del_REFs
    del_df["ALT"] = new_del_ALTs

    snp_df = variants_df[~((variants_df.ALT.str.contains("+", regex=False)) | (variants_df.ALT.str.contains("-", regex=False)))].copy(deep=True).reset_index(drop=True)

    variants_df = pd.concat([snp_df, ins_df, del_df]).sort_values("POS").reset_index(drop=True)

    variants_df.to_csv(output_file, index=False)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-v", "--variantFile", dest='variants_file', type=str, required=True, help='TSV file with VarScan2 variants to be finalized.')
    parser.add_argument("-o", "--outputFile", dest='output_file', type=str, required=True, help='Filepath for output file with finalized variants.')
    
    cmd_line_args = parser.parse_args()
    
    variants_file = cmd_line_args.variants_file
    output_file = cmd_line_args.output_file

    finalize_varscan2(variants_file, output_file)