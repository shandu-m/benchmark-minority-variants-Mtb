#!/home/sm624/.conda/envs/python/bin/python

import pandas as pd

depths = [50, 100, 200, 400, 700]

freqs = [0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]

tools = ["FB", "LF", "MT", "PL", "VD", "VS"]
tool_mapping = {"FB": "FreeBayes", "LF": "LoFreq", "MT": "Mutect2", "PL": "Pilon", "VD": "VarDict", "VS": "VarScan2"}

regions = ["DR", "HT", "LM"]
regions_other = ["DR", "HT", "LM", "other"]
regions_other_nonLM = ["DR", "HT", "other"]

region_mapping = {"DR": "Drug resistance", "HT": "Homopolymer tract", "LM": "Low mappability", "other": "Other"}

region_labels = {"DR": ["DR"], "LM": ["LM", "low_map"], "HT": ["HT"]}
region_labels_lookup = {"DR": "DR", "LM": "LM", "low_map": "LM", "HT": "HT"}

def calc_recall(TP, FN):

    TP, FN = int(TP), int(FN)

    try:
        recall = TP/(TP + FN)
    except ZeroDivisionError:
        recall = 0

    return recall

def calc_precision(TP, FP):

    TP, FP = int(TP), int(FP)

    try:
        precision = TP/(TP + FP)
    except ZeroDivisionError:
        precision = 0

    return precision

def calc_FPR(FP, region_length):

    FP, region_length = int(FP), int(region_length)

    TN = region_length - FP

    try:
        FPR = FP/(FP + TN)
    except ZeroDivisionError:
        FPR = 0

    return FPR

def calc_F1(P, R):

    P, R = float(P), float(R)

    try:
        F1 = (2 * P * R) / (P + R)
    except ZeroDivisionError:
        F1 = 0

    return F1
 
genome_length = 4411532
