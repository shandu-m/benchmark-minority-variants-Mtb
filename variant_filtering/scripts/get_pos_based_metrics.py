import pysam, os
import numpy as np
import pandas as pd
import argparse, vcf, warnings
import subprocess
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sample", dest='sample', type=str, required=True, help='Sample name used to prefix output files.')
parser.add_argument("-b", "--bam", dest='bam_file', type=str, required=True, help='BAM file to get alignment statistics for.')
parser.add_argument("-v", "--variants", dest='variants_file', type=str, required=True, help='CSV file of all variants introduced into a sample/called by at least one tool.')
parser.add_argument("-S", "--simulator", dest='simulator', type=str, required=True, help='Simulator used to simulated sequencing data for this sample.')
parser.add_argument("-G", "--genome", dest='genome', type=str, required=True, help='Reference genome group from which sample was simulated.')

cmd_line_args = parser.parse_args()

sample = cmd_line_args.sample
bam_file = cmd_line_args.bam_file
variants_file = cmd_line_args.variants_file
simulator = cmd_line_args.simulator
genome = cmd_line_args.genome

out_dir = f"/n/scratch/users/s/sm624/FP_characteristics/{simulator}_{genome}/position_metrics"
genome_metrics_dir = f"/n/scratch/users/s/sm624/FP_characteristics/{simulator}_{genome}/genome_metrics"

os.makedirs(out_dir, exist_ok=True)

def compute_mean_base_quality_of_variant_support(bam_file, pos, low_freq_allele):
    
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    base_qualities_low_freq_allele = []

    # default min_base_quality is 13, so bad reads get excluded, which negates the purpose of this function...rolls eyes. min_mapping_quality default is 0
    for pileupcolumn in bam.pileup('NC_000962.3', pos-1, pos, truncate=True, min_base_quality=0):

        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue  # skip deletions and skipped regions

            base = pileupread.alignment.query_sequence[pileupread.query_position]
            qual = pileupread.alignment.query_qualities[pileupread.query_position]

            if base == low_freq_allele:
                base_qualities_low_freq_allele.append(qual)

    return np.mean(base_qualities_low_freq_allele)



############################ GOAL: TO OBTAIN METRICS FOR ALL CALLED VARIANTS ############################


############################ STEP 0: READ IN THE CSV FILE OF ALL VARIANTS ############################

df_variants = pd.read_csv(variants_file)


############################ STEP 1: AVERAGE QUALITY OF BASES SUPPORTING EACH VARIANT ############################

# make empty column first. If df_variants is empty, then Mean_BQ_ALT_allele won't be added, which will cause problems later when trying to access the column
df_variants["MEAN_ALT_BQ"] = np.nan

for i, row in df_variants.iterrows():
    mean_base_qual = compute_mean_base_quality_of_variant_support(bam_file, row['POS'], row['ALT'])
    df_variants.loc[i, "MEAN_ALT_BQ"] = mean_base_qual


############################ STEP 2: OBTAIN STRAND BIAS STATISTICS ############################

SB_tools = ["FB", "MT", "VD", "VS"]

# compute the proportion of variant-supporting reads that are forward sense.
# get the distance from 0.5. Ideally want it to be close to 0.5. So larger the deviation, the less likely the variant is real

for tool in SB_tools:
    df_variants[f"{tool}_STRAND_BIAS"] = np.abs(0.5 - (df_variants[f"{tool}_AO_F"] / df_variants[f"{tool}_AO"]))

############################ STEP 3: ADD COV RATIO, CLIPPED BASES RATIO AND DISCORDANT READS RATIO ############################

genome_metrics_df = pd.read_csv(f"{genome_metrics_dir}/{sample}.genome_metrics.csv").set_index("POS")
genome_metrics_df = genome_metrics_df[["COV_RATIO", "CLIPPED_BASES_RATIO", "DISCORDANT_READS_RATIO"]]

df_variants = df_variants.set_index("POS")

final_variant_df = df_variants.merge(genome_metrics_df, how="left", left_index=True, right_index=True)

final_variant_df.to_csv(f"{out_dir}/{sample}.position_metrics.csv", index=True)



