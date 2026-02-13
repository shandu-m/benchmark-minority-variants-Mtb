#!/home/sm624/bin/anaconda3/bin/python

import pandas as pd
from setup import *
import numpy as np
import sys
import os

def filter_binosnp(tag, unfiltered_vcf_df, output_dir):

    print(">>", tag)

    out_file = open("{0}/{1}.binosnp.filtered.csv".format(output_dir, tag), "w")

    out_file.write("POS,REF,ALT,AF,TYPE\n")
    out_file.flush()
    os.fsync(out_file)

    for i in unfiltered_vcf_df.index:

        # we don't need to do any post-filtering
        # just extract AF

        POS = unfiltered_vcf_df.loc[i, "POS"]
        REF = unfiltered_vcf_df.loc[i, "REF"]
        ALT = unfiltered_vcf_df.loc[i, "ALT"]

        if type(ALT) != float:  
            # handle INDELs
            if "+" in ALT:
                #insertion
                ALT = ALT.strip("+")
                ALT = REF + ALT
            elif "-" in ALT:
                #deletion
                ALT = ALT.strip("-")
                REF = REF + ALT
                ALT = REF[0]
        else:
            ALT = np.nan

        INFO = unfiltered_vcf_df.loc[i, "INFO"].split(";")
        
        TYPE_ix = 0
        TYPE = INFO[TYPE_ix].split("=")[1]
        
        AF_ix = 1
        AF = INFO[AF_ix].split("=")[1]
        if AF == "":
            AF = np.nan
        else:
            AF = float(AF)

        out_file.write("{},{},{},{},{}\n".format(POS, REF, ALT, AF, TYPE))
        out_file.flush()
        os.fsync(out_file)

    out_file.close()

if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("Usage: python filter_binosnp.py [tag] [unfiltered_vcf_path] [output_dir]")
        sys.exit(1)

    tag = sys.argv[1]
    unfiltered_vcf_path = sys.argv[2]
    unfiltered_vcf_df = vcf_to_table(unfiltered_vcf_path)
    output_dir = sys.argv[3]
    
    filter_binosnp(tag, unfiltered_vcf_df, output_dir)
