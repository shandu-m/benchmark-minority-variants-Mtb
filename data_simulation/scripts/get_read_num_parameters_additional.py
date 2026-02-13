#!/home/sm624/bin/anaconda3/bin/python

import pandas as pd
import sys
import numpy as np

def main(tool, depth, freq, var_type):

    read_num_df = pd.read_csv("/home/sm624/projects/mixed_calls/benchmarking/comprehensive/{}/read_num_parameters_additional.csv".format(tool)).set_index(["depth", "freq"])

    read_num = read_num_df.loc[(depth, freq), var_type]
    
    if np.isnan(read_num):
        read_num = ""
    else:
        read_num = int(read_num)

    print(read_num)

if __name__ == "__main__":

    if len(sys.argv) < 5:
        print("Usage: python get_read_num_parameters_additional.py [tool] [depth] [freq] [var_type]")
        sys.exit(1)

    tool = sys.argv[1]
    depth = int(sys.argv[2])
    freq = float(sys.argv[3])
    var_type = sys.argv[4]

    main(tool, depth, freq, var_type)
