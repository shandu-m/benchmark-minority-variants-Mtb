#!/home/sm624/.conda/envs/python/bin/python

import pandas as pd
import numpy as np
import os
import sys
import argparse, warnings

def finalize_freebayes(variants_file):

    variants_df = pd.read_csv(variants_file, sep="\t")

    ### compute AF ###
    variants_df["AF"] = variants_df["AO"].astype(float) / variants_df["DP"].astype(float)
    
    ### label if variant came from a complex type/MNP ###
    complex_type = 1 - variants_df.ORIGIN.isna().astype(int)
    variants_df["FB_COMPLEX_TYPE"] = complex_type

    ### split MNPs (code adapted from Sanjana's function) ###
    df_MNPs = variants_df.query("REF.str.len() == ALT.str.len() & REF.str.len() > 1").copy(deep=True).reset_index(drop=True)
    df_other = variants_df.query("~(REF.str.len() == ALT.str.len() & REF.str.len() > 1)").copy(deep=True).reset_index(drop=True)

    # all the values that must be duplicated
    duplicate_columns = list(set(variants_df.columns) - set(['POS', 'REF', 'ALT']))

    # helps to control ordering so we know that POS, REF, and ALT are at the end
    df_split_MNPs = pd.DataFrame(columns = duplicate_columns + ['POS', 'REF', 'ALT'])
    k = 0

    for _, row in df_MNPs.iterrows():

        split_ref_alleles = list(row['REF'])
        split_alt_alleles = list(row['ALT'])

        # need to update the positions
        start = row['POS']
        pos_lst = np.arange(start, start + len(row['REF']))

        for i, (ref, alt) in enumerate(zip(split_ref_alleles, split_alt_alleles)):

            # some could be the same because the MNP includes both changed and unchanged sites
            if ref != alt:
                df_split_MNPs.loc[k, :] = [row[col] for col in duplicate_columns] + [pos_lst[i], ref, alt]
                k += 1

    variants_df = pd.concat([df_other, df_split_MNPs]).sort_values(['POS']).reset_index(drop=True)
    
    ### determine if variant is supported by multiple haplotypes ###
    hap_counts = variants_df.groupby(["POS", "REF", "ALT"]).size()
    hap_duplicates = hap_counts[hap_counts > 1]
    dup_haps = list(hap_duplicates.index)
    hap_duplicates_df = hap_duplicates.reset_index(drop=False)
    hap_duplicates_df.rename(columns={0:'COUNT'}, inplace=True)
    hap_duplicates_df.set_index("POS", inplace=True)
    
    ### determine if there is evidence for >2 different alleles at a position ###
    pos_counts = variants_df.groupby("POS").size()
    pos_duplicates = pos_counts[pos_counts > 1]
    pos_duplicates_df = pos_duplicates.reset_index(drop=False)
    pos_duplicates_df.rename(columns={0:'COUNT'}, inplace=True)
    pos_duplicates_df.set_index("POS", inplace=True)
    dup_pos = []
    for pos in pos_duplicates_df.index:
        more_occ_than_haps = False
        try:
            more_occ_than_haps = np.all(pos_duplicates_df.loc[pos, "COUNT"] > hap_duplicates_df.loc[pos, "COUNT"])
        except KeyError:
            more_occ_than_haps = True
        if more_occ_than_haps:
            dup_pos.append(pos)
        
    ### label multiple haplotype and multi-alleleic sites ###
    variants_df_indexed = variants_df.set_index(["POS", "REF", "ALT"])
    variants_df_indexed["MULTIPLE_HAPS"] = variants_df_indexed.index.isin(dup_haps).astype(int)
    variants_df = variants_df_indexed.reset_index(drop=False)
    variants_df["MULTIALLELIC"] = variants_df.POS.isin(dup_pos).astype(int)

    origins = variants_df.ORIGIN.values
    variants_df["ORIGIN"] = [int(origin.split(":")[1]) if not pd.isna(origin) else origin for origin in origins]

    GLs = variants_df.GL.values
    ref_GLs = [float(GL.split(",")[0]) if not pd.isna(GL) else GL for GL in GLs]
    alt_GLs = [float(GL.split(",")[1]) if not pd.isna(GL) else GL for GL in GLs]

    variants_df = variants_df.drop("GL", axis=1)
    variants_df["REF_GL"] = ref_GLs
    variants_df["ALT_GL"] = alt_GLs

    variants_df = variants_df.sort_values("POS")

    variants_file_out = variants_file[:-3] + "csv"
    variants_df.to_csv(variants_file_out, index=False)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-v", "--variantFile", dest='variants_file', type=str, required=True, help='TSV file with variants pulled from a normalized FreeBayes VCF to be finalized.')
    
    cmd_line_args = parser.parse_args()
    
    variants_file = cmd_line_args.variants_file

    finalize_freebayes(variants_file)