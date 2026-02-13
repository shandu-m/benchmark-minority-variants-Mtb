#!/home/sm624/.conda/envs/python/bin/python

import os
import io
import pandas as pd

def write(out_file, string):

    out_file.write(string)
    out_file.flush()
    os.fsync(out_file)

def vcf_to_table(path:str, output_type:str='csv', silence:bool=True):

    if output_type in ['csv','c']:

        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]

        if not silence:
            print('<read_vcf> returns: pd.DataFrame')

        return pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                    'QUAL': str, 'FILTER': str, 'INFO': str},
                sep='\t'
                ).rename(columns={'#CHROM': 'CHROM'})

    elif output_type in ['ndarray','n']:

        csv_out = read_vcf(path)

        if not silence:
            print('<read_vcf> returns: ( list, numpy.ndarray(dtype=str) )')
            print('list = column names')
            print('use numpy.column_stack() to reshape as collection of columns')

        return list(csv_out.columns), csv_out.to_numpy(dtype=str)
