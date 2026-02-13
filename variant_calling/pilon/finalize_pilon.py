#!/home/sm624/.conda/envs/python/bin/python

import pandas as pd
import numpy as np
import os
import sys
import argparse, warnings

base_order = {"A":0, "C":1, "G":2, "T":3}
base_order_r = {0:"A", 1:"C", 2:"G", 3:"T"}

def finalize_pilon(variants_file):

    variants_df = pd.read_csv(variants_file, sep="\t")

    ### exclude all records with FILTER tags other than PASS and Amb (confident call) ###
    variants_df = variants_df[variants_df.FILTER.isin(["PASS", "Amb"])]

    ### compute AF ###

    ## SNP records: use QP field to determine AF and base counts ##
    SNP_records_df = variants_df[(variants_df["REF"].str.len() == 1) & (variants_df["ALT"].str.len() <= 1)].copy(deep=True)

    REFs = SNP_records_df.REF.values
    
    QPs = SNP_records_df.QP.values
    AF_arrs = [[float(QP_val)/100 for QP_val in QP.split(",")] for QP in QPs] #use QP as quality-weighted AF

    DPs = SNP_records_df.DP.values
    AC_arrs = [np.array(AF_arr)*DP for AF_arr,DP in zip(AF_arrs, DPs)] #quality-weighted base counts
    
    ALTs = [",".join([base_order_r[b_i] for b_i in base_order_r if AC_arr[b_i] >= 2 and AF_arrs[b_i] >= 0.01 and base_order_r[b_i] != ref]) for ref, AC_arr, AF_arrs in zip(REFs, AC_arrs, AF_arrs)]
    ACs = [",".join([str(AC_arr[b_i]) for b_i in base_order_r if AC_arr[b_i] >= 2 and AF_arrs[b_i] >= 0.01 and base_order_r[b_i] != ref]) for ref, AC_arr, AF_arrs in zip(REFs, AC_arrs, AF_arrs)]
    new_AFs = [",".join([str(AF_arrs[b_i]) for b_i in base_order_r if AC_arr[b_i] >= 2 and AF_arrs[b_i] >= 0.01 and base_order_r[b_i] != ref]) for ref, AC_arr, AF_arrs in zip(REFs, AC_arrs, AF_arrs)]

    SNP_records_df["ALT"] = ALTs
    SNP_records_df["AC"] = ACs
    SNP_records_df["AF"] = new_AFs

    # keep positions with a valid ALT allele
    SNP_records_df = SNP_records_df[SNP_records_df.ALT.str.len() > 0].reset_index(drop=True)

    # split multi-allelic positions
    nonMA_SNP_records_df = SNP_records_df[~SNP_records_df.ALT.str.contains(",")].copy(deep=True).reset_index(drop=True)

    MA_SNP_records_df = SNP_records_df[SNP_records_df.ALT.str.contains(",")].copy(deep=True).reset_index(drop=True)
    
    df_add_i = MA_SNP_records_df.shape[0]

    for i in MA_SNP_records_df.index:
    
        alts = MA_SNP_records_df.loc[i, "ALT"].split(",")
        acs = MA_SNP_records_df.loc[i, "AC"].split(",")
        afs = MA_SNP_records_df.loc[i, "AF"].split(",")
        
        MA_SNP_records_df.loc[i, "ALT"] = alts[0]
        MA_SNP_records_df.loc[i, "AC"] = acs[0]
        MA_SNP_records_df.loc[i, "AF"] = afs[0]
        
        for alt_i in range(1, len(alts)):
            MA_SNP_records_df.loc[df_add_i] = MA_SNP_records_df.loc[i]
            MA_SNP_records_df.loc[df_add_i, "ALT"] = alts[alt_i]
            MA_SNP_records_df.loc[df_add_i, "AC"] = acs[alt_i]
            MA_SNP_records_df.loc[df_add_i, "AF"] = afs[alt_i]
            df_add_i += 1

    SNP_records_df = pd.concat([nonMA_SNP_records_df, MA_SNP_records_df])
    SNP_records_df["AC"] = SNP_records_df["AC"].astype(float)
    SNP_records_df["AF"] = SNP_records_df["AF"].astype(float)

    ## non-SNP records: not all have DP and AF populated by Pilon ##
    SV_records_df = variants_df[~((variants_df["REF"].str.len() == 1) & (variants_df["ALT"].str.len() <= 1))].copy(deep=True)

    # keep SV records with no DP or AF reported as is
    NaN_SV_records_df = SV_records_df[SV_records_df.AF.isna()].copy(deep=True)

    # filter AC and AF for SVs with these fields
    nonNaN_SV_records_df = SV_records_df[~SV_records_df.AF.isna()].copy(deep=True)

    DPs = nonNaN_SV_records_df.DP.values
    AFs = nonNaN_SV_records_df.AF.values

    ACs = [AF*DP for AF,DP in zip(AFs, DPs)]
    nonNaN_SV_records_df["AC"] = ACs

    filtered_nonNaN_SV_records_df = nonNaN_SV_records_df[(nonNaN_SV_records_df.AC >= 2) & (nonNaN_SV_records_df.AF >= 0.01)]

    SV_records_df = pd.concat([NaN_SV_records_df, filtered_nonNaN_SV_records_df])

    variants_df = pd.concat([SNP_records_df, SV_records_df])

    variants_df = variants_df.sort_values("POS")
    
    variants_file_out = variants_file[:-3] + "csv"
    variants_df.to_csv(variants_file_out, index=False)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-v", "--variantFile", dest='variants_file', type=str, required=True, help='TSV file with variants pulled from a normalized Pilon VCF to be finalized.')
    
    cmd_line_args = parser.parse_args()
    
    variants_file = cmd_line_args.variants_file

    finalize_pilon(variants_file)