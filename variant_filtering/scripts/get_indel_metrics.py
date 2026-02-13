import numpy as np
import pandas as pd
import os, glob, warnings, pysam
warnings.filterwarnings('ignore')
from Bio import Seq, SeqIO
import scipy.stats as st
import argparse

h37Rv_path = "/n/data1/hms/dbmi/farhat/Sanjana/H37Rv"
h37Rv_seq = SeqIO.read(os.path.join(h37Rv_path, "GCF_000195955.2_ASM19595v2_genomic.gbff"), "genbank")
h37Rv_genes = pd.read_csv(os.path.join(h37Rv_path, "mycobrowser_h37rv_genes_v4.csv"))
h37Rv_regions = pd.read_csv(os.path.join(h37Rv_path, "mycobrowser_h37rv_v4.csv"))

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sample", dest='sample', type=str, required=True, help='Sample name used to prefix output files.')
parser.add_argument("-b", "--bam", dest='bam_file', type=str, required=True, help='BAM file to get alignment statistics for.')
parser.add_argument("-v", "--indelVariants", dest='variants_file', type=str, required=True, help='CSV file of all candidate INDEL variants.')
parser.add_argument("-o", "--outDir", dest='out_dir', type=str, required=True, help='Output directory in which to put INDELs with adjusted AFs.')
parser.add_argument("-S", "--simulator", dest='simulator', type=str, required=True, help='Simulator used to simulated sequencing data for this sample.')
parser.add_argument("-G", "--genome", dest='genome', type=str, required=True, help='Reference genome group from which sample was simulated.')
parser.add_argument("-c", "--chromosome", dest='chromosome', type=str, required=True, help='Chromosome name used for alignment.')

cmd_line_args = parser.parse_args()

sample = cmd_line_args.sample
bam_file = cmd_line_args.bam_file
variants_file = cmd_line_args.variants_file
out_dir = cmd_line_args.out_dir
simulator = cmd_line_args.simulator
genome = cmd_line_args.genome
chromosome = cmd_line_args.chromosome

def compute_adjusted_AF_indels(sample, df_indels, chrom, bam_file, n_prop=0.50, dist_from_indel=10):
    
    df_indels = df_indels.reset_index(drop=True)
    
    bam = pysam.AlignmentFile(bam_file, "rb")

    for i, row in df_indels.iterrows():

        pos = row['POS']
        ref = row['REF']
        alt = row['ALT']

        indel_length = abs(len(ref) - len(alt))

        if len(alt) > len(ref):
            search_string = 'I'
            search_val = 1
        else:
            search_string = 'D'
            search_val = 2

        reads_overlapping = []

        for read in bam.fetch(chrom, pos - 1, pos):  # pysam uses 0-based, half-open

            # exclude reads with supplementary alignments, secondary alignments, or are discordantly paired
            if not read.is_supplementary and not read.is_secondary and read.is_proper_pair:
               
                if read.cigarstring is not None:
                    # compute the total number of I or D values in the cigar string. So if the cigar is 1M10I4X15I, and search_string = 'I', it returns 25
                    # if there are multiple indels at the same site, they will be reflected by different values of the num_indel column and will split the support accordingly
                    indel_length_from_cigar = sum(length for (op, length) in (read.cigartuples or []) if op == search_val and length > 0)
                        
                    # get number of bases soft clipped from this read
                    num_soft_clipped_bases = sum(length for (op, length) in (read.cigartuples or []) if op == 4)

                    reads_overlapping.append({
                        "read_name": read.query_name,
                        "ref_start": read.reference_start + 1,  # make it 1-based
                        "ref_end": read.reference_end,
                        "num_indel": indel_length_from_cigar,
                        "num_soft_clipped": num_soft_clipped_bases,
                        "cigar": read.cigarstring,
                        "cigartuples": read.cigartuples,
                    })

        if len(reads_overlapping) > 0:

            df_overlapping_reads = pd.DataFrame(reads_overlapping).sort_values(['ref_start', 'ref_end'])

            # take the middle N% of reads? No, this is hard to tune because it needs to be low enough to work for some fixed indels, 
            # but high enough for very rare indels (like AF < 10%) so that we still detect them
            # df_reads_close_to_indel = keep_middle_percentage_of_reads(df_overlapping_reads, n_prop)

            # include a read only if neither its start nor end are within 10 base pairs of the indel site
            df_reads_close_to_indel = df_overlapping_reads.loc[(abs(df_overlapping_reads['ref_start'] - pos) >= dist_from_indel) & 
                                                               (abs(df_overlapping_reads['ref_end'] - pos) >= dist_from_indel)
                                                              ]

            # also exclude reads that start or end in the middle of the indel            
            # this position is inclusive. So the read must not have a start or end within this
            indel_end = pos + indel_length

            df_reads_close_to_indel = df_reads_close_to_indel.loc[~((df_reads_close_to_indel['ref_start'] >= pos) & (df_reads_close_to_indel['ref_start'] <= indel_end))]

            df_reads_close_to_indel = df_reads_close_to_indel.loc[~((df_reads_close_to_indel['ref_end'] >= pos) & (df_reads_close_to_indel['ref_end'] <= indel_end))]

            # AND include a read only if neither its start nor end are within 10 base pairs of the deletion END
            df_reads_close_to_indel = df_reads_close_to_indel.loc[(abs(df_reads_close_to_indel['ref_start'] - indel_end) >= dist_from_indel) & 
                                                                  (abs(df_reads_close_to_indel['ref_end'] - indel_end) >= dist_from_indel)
                                                                 ]


            # AND exclude reads with ANY soft clipping from contributing to the depth
            df_reads_close_to_indel = df_reads_close_to_indel.query("num_soft_clipped == 0")
            
            # these are very low quality indels. If all reads have been excluded, then it was probably an area of tons of discordant alignments, meaning there's no low frequency indel. 
            # It's reference bias due to SVs
            if len(df_reads_close_to_indel) == 0:
                df_indels.loc[i, ['Indel_Support', 'Total_Reads', 'AF_Adj']] = [0, 0, 0]
            else:
                num_reads_close_supporting_indel = len(df_reads_close_to_indel.query("num_indel == @indel_length"))

                adj_AF = num_reads_close_supporting_indel / len(df_reads_close_to_indel)

                # first: there are reads supporting the indel (num_indel.max() > 0) and reads not supporting it (num_indel.min() == 0)
                # need to make separate cases for if there are multiple indels at the site, which happens in repetitive regions
                # where there may be one or multiple copies of the repeat unit inserted or deleted
                # also need to check if indel_length is in df_reads_close_to_indel.num_indel.values because the reads supporting it may be too close to the indel and therefore removed above
                if df_reads_close_to_indel.num_indel.min() == 0 and df_reads_close_to_indel.num_indel.max() > 0 and indel_length in df_reads_close_to_indel.num_indel.values:
                   
                    # add a check to see if the start and end of reads that support and don't support the indel are significantly different. If they are, then it's probably a fixed indel.
                    # KS test is very sensitive, using the median (Mann-Whitney U test) is more robust to outliers
                    start_pval = st.mannwhitneyu(df_reads_close_to_indel.query("num_indel==@indel_length").ref_start,
                                                df_reads_close_to_indel.query("num_indel==0").ref_start
                                               ).pvalue

                    end_pval = st.mannwhitneyu(df_reads_close_to_indel.query("num_indel==@indel_length").ref_end,
                                                df_reads_close_to_indel.query("num_indel==0").ref_end
                                               ).pvalue

                    # set the AF as 1
                    if start_pval < 0.01 and end_pval < 0.01:

                        # the reads that don't support the indel just don't fully span the indel. They look like they don't support the indel but it's artifactual
                        # so remove them and recalculate the AF. These reads shouldn't contribute to the depth
                        df_reads_close_to_indel = df_reads_close_to_indel.query("num_indel > 0")

                        #  These will not support the indel, which doesn't mean that the indel is present at low frequency. It's artifactual
                        df_indels.loc[i, ['Indel_Support', 'Total_Reads', 'AF_Adj']] = [num_reads_close_supporting_indel, 
                                                                                        len(df_reads_close_to_indel),
                                                                                        num_reads_close_supporting_indel / len(df_reads_close_to_indel)
                                                                                       ]

                    else:
                        df_indels.loc[i, ['Indel_Support', 'Total_Reads', 'AF_Adj']] = [num_reads_close_supporting_indel, 
                                                                                        len(df_reads_close_to_indel),
                                                                                        adj_AF
                                                                                       ]

                # there are only reads supporting the indel. Then no need for the KS test
                else:
                    df_indels.loc[i, ['Indel_Support', 'Total_Reads', 'AF_Adj']] = [num_reads_close_supporting_indel, 
                                                                                    len(df_reads_close_to_indel),
                                                                                    adj_AF
                                                                                   ]

        else:
            # there should always be overlapping reads. If not, it means there's an SV or something with few to no reads in the pileup
            # but if that's the case, then there is def no indel, so remove these from the unfixed indels table
            df_indels.loc[i, ['Indel_Support', 'Total_Reads', 'AF_Adj']] = [np.nan, np.nan, np.nan]

    bam.close()

    return df_indels

############################ GOAL: TO COMPUTE AN ADJUSTED AF FOR INDELS ############################

dist_from_indel = 10
        
df_indels = pd.read_csv(f"{variants_file}")

df_indels = df_indels.sort_values('POS')

if len(df_indels) > 0:
    
    df_indels_adjusted_AF = compute_adjusted_AF_indels(sample,
                                                       df_indels, 
                                                       chromosome, 
                                                       f"{bam_file}", 
                                                       dist_from_indel=dist_from_indel,
                                                      )
    
else:

    df_indels_adjusted_AF = pd.DataFrame(columns=['POS', 'REF', 'ALT', 'GT', 'REGION', 'FB_AF', 'AO', 'FB_AO_F',
                                                  'FB_AO_R', 'FB_DP', 'FB_MQ', 'FB_QUAL', 'FB_ORIGIN', 'FB_COMPLEX_TYPE',
                                                  'FB_MULTIPLE_HAPS', 'FB_MULTIALLELIC', 'FB_REF_GL', 'FB_ALT_GL',
                                                  'Indel_Support', 'Total_Reads', 'AF_Adj'])
    
df_indels_adjusted_AF.to_csv(f"{out_dir}/{sample}.INDELs.csv", index=False)










