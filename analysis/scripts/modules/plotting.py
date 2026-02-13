#!/home/sm624/.conda/envs/python/bin/python

import matplotlib.pyplot as plt

LF_col, VS_col, VD_col, MT_col, FB_col, PL_col = "#ece133", "#d55e00", "#949494", "#fbafe4", "#56b4e9", "#029e73"
tool_colors = {"LF": LF_col, "VS": VS_col, "VD": VD_col, "MT": MT_col, "FB": FB_col, "PL": PL_col}
tool_colors_ = {"LoFreq": LF_col, "VarScan2": VS_col, "VarDict": VD_col, "Mutect2": MT_col, "FreeBayes": FB_col, "Pilon": PL_col}

depth_pal = ["#1e90ff", "#2065ab", "#215080", "#223b56", "#23262c"]
depth_colors = {50: depth_pal[0], 100: depth_pal[1], 200: depth_pal[2], 400: depth_pal[3], 700: depth_pal[4]}

#DR_col, LM_col, HT_col, other_col = "#cc78bc", "red", "#de8f05", "grey"
DR_col, LM_col, HT_col, other_col = "#0000a2", "#bc272d", "#e9c716", "grey"
region_colors = {"DR": DR_col, "HT": HT_col, "LM": LM_col, "other": other_col}
region_colors_ = {"Drug resistance": DR_col, "Homopolymer tract": HT_col, "Low mappability": LM_col, "Other": other_col}

Rv_L1234_palette = {"H37Rv": "#1E90FF", "L1234": "#FFD700"}
Rv_L1234_palette_ = {"H37Rv": "#1E90FF", "L1-4": "#FFD700"}

TP_col, FP_col = "#1a80bb", "#ea801c"
P_palette = {"TP": TP_col, "FP": FP_col}