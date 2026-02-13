#!/home/sm624/.conda/envs/python/bin/python

import pandas as pd
import pickle
import numpy as np

# drug resistance
DR_pos_geneIDName_pkl = "/home/sm624/projects/mixed_calls/benchmarking/comprehensive/analysis/all/FINAL/data/DR_pos_geneIDName_dict.pkl"
with open(DR_pos_geneIDName_pkl, "rb") as in_f:
    DR_pos = pickle.load(in_f)

# homopolymer tract
HT_pos_varBase_dict = "/home/sm624/projects/mixed_calls/benchmarking/comprehensive/analysis/all/FINAL/data/HT_pos_varBase_dict.pkl"
with open(HT_pos_varBase_dict, "rb") as in_f:
    HT_pos = pickle.load(in_f)

# comprehensive LM regions
comprehensive_LM_regions = "/n/data1/hms/dbmi/farhat/shandu/projects/TRUST/data/L1234.H37Rv.pupmapper.indelHomologRegions.RLC.geneHomologs.first500bp.bed"

LM_regions_df = pd.read_csv(comprehensive_LM_regions, names=["chrom", "start", "stop"], sep="\t")
LM_pos = []
for i in LM_regions_df.index:
    start, stop = LM_regions_df.loc[i, "start"], LM_regions_df.loc[i, "stop"]
    LM_pos += list(np.arange(start+1, stop+1))

# simulation LM regions
LM_pos_list_pkl = "/home/sm624/projects/mixed_calls/benchmarking/comprehensive/analysis/all/FINAL/data/LM_pos_list.pkl"
with open(LM_pos_list_pkl, "rb") as in_f:
    sim_LM_pos = pickle.load(in_f)

# function to determine position type
def get_pos_type(pos, DR_pos=DR_pos, HT_pos=HT_pos, LM_pos=LM_pos):

    if pos in DR_pos:
        return "DR"
    elif pos in HT_pos:
        return "HT"
    elif pos in LM_pos:
        return "LM"
    return "other"

DR_region_length = len(DR_pos)
HT_region_length = len(HT_pos)
LM_region_length = len(LM_pos)
sim_LM_region_length = len(sim_LM_pos)
