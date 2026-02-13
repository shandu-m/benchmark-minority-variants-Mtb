#!/home/sm624/.conda/envs/vcflib/bin/python

import pandas as pd
import numpy as np
import sys
import pickle
import argparse

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

def get_pos_type_sim(pos, DR_pos=DR_pos, HT_pos=HT_pos, LM_pos=sim_LM_pos):

    if pos in DR_pos:
        return "DR"
    elif pos in HT_pos:
        return "HT"
    elif pos in LM_pos:
        return "LM"
    return "other"

def make_summary(tag, genome, source_dir, output_dir):

    ## ground truth mutations ##
    if genome == "L1234":

        GT_df = pd.read_csv("/home/sm624/projects/mixed_calls/benchmarking/comprehensive/lineage123/mutant_summaries_final/L123_mutant.csv")
        
        GT_df["GT"] = 1
        GT_df_indexed = GT_df.set_index(["POS", "REF", "ALT"]).drop("TYPE", axis=1)
    
        # if genome == L1234 then add ground truth lineage mutations too
        lineage = tag.split("_")[1]
        baseline_df = pd.read_csv(f"/home/sm624/projects/mixed_calls/benchmarking/comprehensive/lineage123/baseline_variants/{lineage}_baseline.csv")
        
        baseline_df["L_GT"] = 1
        baseline_df_indexed = baseline_df.set_index(["POS", "REF", "ALT"])
        
        GT_df_indexed = GT_df_indexed.merge(baseline_df_indexed, how="outer", left_index=True, right_index=True)

    else:

        # H37Rv
        mutant_num = tag.split("_")[3]
        GT_df = pd.read_csv(f"/home/sm624/projects/mixed_calls/benchmarking/comprehensive/mutant_summaries_final/mutant_{mutant_num}.csv")
        
        GT_df["GT"] = 1
        GT_df_indexed = GT_df.set_index(["POS", "REF", "ALT"]).drop("TYPE", axis=1)

    ## FreeBayes variants ##
    FB_df = pd.read_csv(f"{source_dir}/{tag}.freebayes.csv")

    # handle multiple occurences from different haplotypes
    FB_indexed_ = FB_df.set_index(["POS", "REF", "ALT"])
    FB_dup_variants = FB_indexed_[FB_indexed_.index.duplicated()].index
    FB_ix_to_drop = []
    for dup_var in FB_dup_variants:
        dup_entry = FB_df[FB_indexed_.index == dup_var]
        dup_entry_ix = dup_entry.index
        keep_ix = dup_entry["AF"].idxmax()
        drop_ix = list(dup_entry_ix.drop(keep_ix))
        FB_ix_to_drop += drop_ix
    FB_updated_df = FB_df.drop(FB_ix_to_drop)

    FB_updated_df["FB_found"] = 1
    FB_updated_df = FB_updated_df.rename({"AF":"FB_AF", "AO":"FB_AO", "AO_F":"FB_AO_F", "AO_R":"FB_AO_R",
                                          "DP":"FB_DP", "QUAL":"FB_QUAL", "MQM":"FB_MQ", "ORIGIN":"FB_ORIGIN",
                                          "MULTIPLE_HAPS":"FB_MULTIPLE_HAPS", "MULTIALLELIC":"FB_MULTIALLELIC",
                                          "REF_GL":"FB_REF_GL", "ALT_GL":"FB_ALT_GL"}, axis=1).drop("ID", axis=1)
    FB_df_indexed = FB_updated_df.set_index(["POS", "REF", "ALT"])

    ## LoFreq variants ##
    LF_df = pd.read_csv(f"{source_dir}/{tag}.lofreq.csv")

    LF_df["LF_found"] = 1
    LF_df = LF_df.rename({"QUAL":"LF_QUAL", "DP":"LF_DP", "AF":"LF_AF", "AC":"LF_AO"}, axis=1)
    LF_df_indexed = LF_df.set_index(["POS", "REF", "ALT"])

    ## Mutect2 variants ##
    MT_df = pd.read_csv(f"{source_dir}/{tag}.mutect2.csv")

    # handle multiple occurences from different haplotypes
    MT_indexed_ = MT_df.set_index(["POS", "REF", "ALT"])
    MT_dup_variants = MT_indexed_[MT_indexed_.index.duplicated()].index
    MT_ix_to_drop = []
    for dup_var in MT_dup_variants:
        dup_entry = MT_df[MT_indexed_.index == dup_var]
        dup_entry_ix = dup_entry.index
        keep_ix = dup_entry["AF"].idxmax()
        drop_ix = list(dup_entry_ix.drop(keep_ix))
        MT_ix_to_drop += drop_ix
    MT_updated_df = MT_df.drop(MT_ix_to_drop)

    MT_updated_df["MT_found"] = 1

    MT_updated_df["MQ"] = [int(x.split(",")[1]) for x in MT_updated_df.MQ.values]
    
    AO_F = [int(x.split("|")[1].split(",")[0]) for x in MT_updated_df.SB.values]
    AO_R = [int(x.split("|")[1].split(",")[1]) for x in MT_updated_df.SB.values]
    
    MT_updated_df["MT_AO_F"] = AO_F
    MT_updated_df["MT_AO_R"] = AO_R
    
    MT_updated_df = MT_updated_df.rename({"AF":"MT_AF", "FILTER":"MT_FILTER", "DP":"MT_DP", "MQ":"MT_MQ",
                                          "ORIGIN":"MT_ORIGIN", "PID":"MT_PID", "AC":"MT_AO",
                                          "MULTIPLE_HAPS":"MT_MULTIPLE_HAPS",
                                          "MULTIALLELIC":"MT_MULTIALLELIC"}, axis=1).drop(["ID", "SB"], axis=1)
    MT_df_indexed = MT_updated_df.set_index(["POS", "REF", "ALT"])
    
    ## Pilon variants ##
    PL_df = pd.read_csv(f"{source_dir}/{tag}.pilon.csv")

    PL_df["PL_found"] = 1
    PL_df = PL_df.rename({"QUAL":"PL_QUAL", "FILTER":"PL_FILTER", "DP":"PL_DP", "MQ":"PL_MQ", 
                          "AF":"PL_AF", "AC":"PL_AO"}, axis=1).drop("QP", axis=1) 
    PL_df_indexed = PL_df.set_index(["POS", "REF", "ALT"])
    
    ## VarDict variants ##
    VD_df = pd.read_csv(f"{source_dir}/{tag}.vardict.csv")

    # handle multiple occurences from different haplotypes
    VD_indexed_ = VD_df.set_index(["POS", "REF", "ALT"])
    VD_dup_variants = VD_indexed_[VD_indexed_.index.duplicated()].index
    VD_ix_to_drop = []
    for dup_var in VD_dup_variants:
        dup_entry = VD_df[VD_indexed_.index == dup_var]
        dup_entry_ix = dup_entry.index
        keep_ix = dup_entry["AF"].idxmax()
        drop_ix = list(dup_entry_ix.drop(keep_ix))
        VD_ix_to_drop += drop_ix
    VD_updated_df = VD_df.drop(VD_ix_to_drop)

    VD_updated_df["VD_found"] = 1

    AO_F = [int(x.split(":")[0]) for x in VD_updated_df.VARBIAS.values]
    AO_R = [int(x.split(":")[1]) for x in VD_updated_df.VARBIAS.values]
    
    VD_updated_df["VD_AO_F"] = AO_F
    VD_updated_df["VD_AO_R"] = AO_R
    
    VD_updated_df = VD_updated_df.rename({"FILTER":"VD_FILTER", "AF":"VD_AF", "DP":"VD_DP", "MQ":"VD_MQ", 
                                          "ORIGIN":"VD_ORIGIN", "AC":"VD_AO", "MULTIPLE_HAPS":"VD_MULTIPLE_HAPS",
                                          "MULTIALLELIC":"VD_MULTIALLELIC"}, axis=1).drop(["REFBIAS", "VARBIAS"], axis=1)
    VD_df_indexed = VD_updated_df.set_index(["POS", "REF", "ALT"])

    ## VarScan2 variants ##
    VS_df = pd.read_csv(f"{source_dir}/{tag}.varscan2.csv")
    
    VS_df["VS_found"] = 1

    VS_df = VS_df.rename({"AF":"VS_AF", "DP":"VS_DP", "AO_F":"VS_AO_F", "AO_R":"VS_AO_R",
                          "AC":"VS_AO"}, axis=1)
    VS_df_indexed = VS_df.set_index(["POS", "REF", "ALT"])


    ## MERGE ##
    variant_summary_df = GT_df_indexed.merge(FB_df_indexed, how="outer", left_index=True, right_index=True).merge(LF_df_indexed, how="outer", left_index=True, right_index=True).merge(MT_df_indexed, how="outer", left_index=True, right_index=True).merge(PL_df_indexed, how="outer", left_index=True, right_index=True).merge(VD_df_indexed, how="outer", left_index=True, right_index=True).merge(VS_df_indexed, how="outer", left_index=True, right_index=True)

    variant_summary_df = variant_summary_df.reset_index(drop=False).sort_values("POS")

    # fill NaN values in with 0 for tool_found and GT columns and L_GT if L1234
    tools = ["FB", "LF", "MT", "PL", "VD", "VS"]
    for tool in tools:
        variant_summary_df[f"{tool}_found"] = variant_summary_df[f"{tool}_found"].fillna(0)

    variant_summary_df["GT"] = variant_summary_df["GT"].fillna(0)

    if genome == "L1234":
        variant_summary_df["L_GT"] = variant_summary_df["L_GT"].fillna(0)

    # annotate with position
    variant_summary_df["REGION"] = variant_summary_df["POS"].apply(get_pos_type)
    variant_summary_df["REGION_sim"] = variant_summary_df["POS"].apply(get_pos_type_sim)

    # fix types
    variant_summary_df["POS"] = variant_summary_df["POS"].astype(int)
    variant_summary_df["REF"] = variant_summary_df["REF"].astype(str)
    variant_summary_df["ALT"] = variant_summary_df["ALT"].astype(str)
    
    variant_summary_df["REGION"] = variant_summary_df["REGION"].astype(str)
    variant_summary_df["REGION_sim"] = variant_summary_df["REGION_sim"].astype(str)
    variant_summary_df["GT"] = variant_summary_df["GT"].astype(int)
      
    for tool in tools:
        variant_summary_df["{}_found".format(tool)] = variant_summary_df["{}_found".format(tool)].astype(int)
        variant_summary_df["{}_AF".format(tool)] = variant_summary_df["{}_AF".format(tool)].astype(float)

    if genome == "L1234":
        variant_summary_df["L_GT"] = variant_summary_df["L_GT"].astype(int)

    variant_summary_df.to_csv(f"{output_dir}/{tag}.variant_summary.csv", index=False)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-t", "--tag", dest='tag', type=str, required=True, help='Sample name.')
    parser.add_argument("-g", "--genome", dest='genome', type=str, required=True, help='Genome from which sample was simulated (H37Rv/L1234).')
    parser.add_argument("-s", "--source_dir", dest='source_dir', type=str, required=True, help='Directory in which to find finalized variant files for each tool.')
    parser.add_argument("-o", "--output_dir", dest='output_dir', type=str, required=True, help='Directory in which to write variant summary file for sample.')
    
    cmd_line_args = parser.parse_args()

    tag = cmd_line_args.tag
    genome = cmd_line_args.genome
    source_dir = cmd_line_args.source_dir
    output_dir = cmd_line_args.output_dir

    make_summary(tag, genome, source_dir, output_dir)