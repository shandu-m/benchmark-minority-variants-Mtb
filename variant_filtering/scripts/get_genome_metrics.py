import pysam, os
import numpy as np
import pandas as pd
import argparse, vcf, warnings
import subprocess
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sample", dest='sample', type=str, required=True, help='Sample name used to prefix output files.')
parser.add_argument("-b", "--bam", dest='bam_file', type=str, required=True, help='BAM file to get alignment statistics for.')
parser.add_argument("-S", "--simulator", dest='simulator', type=str, required=True, help='Simulator used to simulated sequencing data for this sample.')
parser.add_argument("-G", "--genome", dest='genome', type=str, required=True, help='Reference genome group from which sample was simulated.')
parser.add_argument("-g", "--genomeFile", dest='genome_file', type=str, default="/home/sak0914/MtbLongitudinalDiversity/lowAF_variant_calling/references/ref_genome/H37Rv_NC_000962.3.fna", help='Reference genome FASTA file')

cmd_line_args = parser.parse_args()

sample = cmd_line_args.sample
bam_file = cmd_line_args.bam_file
simulator = cmd_line_args.simulator
genome = cmd_line_args.genome
genome_file = cmd_line_args.genome_file

out_dir = f"/n/scratch/users/s/sm624/FP_characteristics/{simulator}_{genome}/genome_metrics"

os.makedirs(out_dir, exist_ok=True)

soft_clipping_BED_file = os.path.join(out_dir, f"{sample}.softClips.bed")
soft_clips_by_pos_file = os.path.join(out_dir, f"{sample}.softClips.tsv.gz")

discordant_alns_BED_file = os.path.join(out_dir, f"{sample}.DiscordantReads.bed")
discordant_alns_by_pos_file = os.path.join(out_dir, f"{sample}.DiscordantReads.tsv.gz")

h37Rv_path = "/n/data1/hms/dbmi/farhat/Sanjana/H37Rv"
h37Rv_regions = pd.read_csv(os.path.join(h37Rv_path, "mycobrowser_h37rv_v4.csv"))

def fasta_length(genome_file):

    genome_length = 0
    
    with open(genome_file, 'r') as in_f:

        header = in_f.readline()
        
        line = in_f.readline().strip()
        while line:
            genome_length += len(line)
            line = in_f.readline().strip()

    return genome_length


def get_orientation(read):
    """Return orientation of the pair: LR, RL, LL, RR, or None if not valid pair."""
    if not (read.is_paired and not read.is_unmapped and not read.mate_is_unmapped):
        return None

    # different chromosomes
    if read.reference_id != read.next_reference_id:
        return "DIFF_CHR"

    mate_start = read.next_reference_start
    this_start = read.reference_start

    # Determine left and right mate
    if this_start <= mate_start:
        left_rev  = read.is_reverse
        right_rev = read.mate_is_reverse
    else:
        left_rev  = read.mate_is_reverse
        right_rev = read.is_reverse

    if not left_rev and right_rev:
        return "LR"   # forward (left) + reverse (right) = normal Illumina
    if left_rev and not right_rev:
        return "RL"   # reverse (left) + forward (right) = IGV green
    if not left_rev and not right_rev:
        return "LL"   # both forward
    if left_rev and right_rev:
        return "RR"   # both reverse


def save_table_of_soft_clips(sample, bam_file, ref_genome_file):
        
    genome_length = fasta_length(ref_genome_file)
    
    bam = pysam.AlignmentFile(bam_file, "rb")
    chrom = bam.references[0]

    records = []
    for read in bam.fetch():

        # skip unmapped reads
        if read.is_unmapped:
            continue
            
        num_left_clip = 0
        num_right_clip = 0

        if read.cigartuples:
            
            # check first operator
            if read.cigartuples[0][0] == 4:  # 4 = soft clip
                num_left_clip = read.cigartuples[0][1]

            # check last operator
            if read.cigartuples[-1][0] == 4:
                num_right_clip = read.cigartuples[-1][1]
                
        # add 1 because it's 0-indexed half-open
        start = read.reference_start + 1

        # don't add 1 because it's half-open at the end
        end = read.reference_end

        records.append({
            "read_name": read.query_name,
            "chrom": chrom,
            "start": start,
            "end": end,
            "cigar_string": read.cigarstring,
            "soft_clipped_bases": num_left_clip + num_right_clip, # total
            "left_clip_start": start - num_left_clip if num_left_clip > 0 else np.nan,
            "left_clip_end": start - 1 if num_left_clip > 0 else np.nan,
            "right_clip_start": end + 1 if num_right_clip > 0 else np.nan,
            "right_clip_end": end + num_right_clip if num_right_clip > 0 else np.nan,
        })

    bam.close()

    df = pd.DataFrame(records)
    
    # there will be clipping at the ends because the genome is circular, but we aligned reads to the linearized genome
    # exclude those from this computation because they're more artifactual
    # THESE ARE 1-INDEXED CLOSED INTERVALS, so the genome ranges from 1 - genome_length, INCLUSIVE
    # so exclude any soft-clipped regions that extend outside of the ref genome AT ALL
    df_soft_clips = df.dropna(subset=['left_clip_start', 'right_clip_start'], how='all').query("~(left_clip_start < 1) & ~(right_clip_end > @genome_length)").reset_index(drop=True)
    
    # convert this to a BED file, then use bedtools genomecov because we're getting a coverage track of soft clipping
    # make sure to get both left and right clips. Expand them like this
    df_soft_clips_BED = pd.concat([df_soft_clips.dropna(subset='left_clip_start')[['chrom', 'left_clip_start', 'left_clip_end']].rename(columns={'left_clip_start': 'BEG', 'left_clip_end': 'END'}),
                                   df_soft_clips.dropna(subset='right_clip_start')[['chrom', 'right_clip_start', 'right_clip_end']].rename(columns={'right_clip_start': 'BEG', 'right_clip_end': 'END'})
                                  ])

    # because we expanded, the length should be at least as long as df_soft_clips. And more than likely longer because some reads have both left and right soft clipping
    assert len(df_soft_clips_BED) >= len(df_soft_clips)
    # print(len(df_soft_clips_BED), len(df_soft_clips))

    # convert to 0-indexed half open intervals
    df_soft_clips_BED['BEG'] -= 1
    df_soft_clips_BED[['BEG', 'END']] = df_soft_clips_BED[['BEG', 'END']].astype(int)
    
    df_soft_clips_BED.to_csv(soft_clipping_BED_file, sep='\t', header=None, index=False)
    
    ref_genome_chrom_lengths_file = f"{ref_genome_file.replace('.fna', '')}.chrom_lengths.txt"
    
    # then run the following to compute the number of soft clipped 
    if not os.path.isfile(ref_genome_chrom_lengths_file):
        subprocess.run(f"cut -f1,2 {ref_genome_file}.fai > {ref_genome_chrom_lengths_file}", shell=True)

    subprocess.run(f"bedtools genomecov -i {soft_clipping_BED_file} -g {ref_genome_chrom_lengths_file} -d | gzip -c > {soft_clips_by_pos_file}", shell=True)
    
    


def save_table_of_discordant_reads(sample, bam_file, ref_genome_file):
        
    genome_length = fasta_length(ref_genome_file)
    
    bam = pysam.AlignmentFile(bam_file, "rb")
    chrom = bam.references[0]

    records = []

    for read in bam.fetch():

        # skip unmapped reads
        if read.is_unmapped:
            continue

        orientation = get_orientation(read)

        if orientation != 'LR':
            
            # add 1 because it's 0-indexed half-open
            start = read.reference_start + 1

            # don't add 1 because it's half-open at the end
            end = read.reference_end
        
            records.append({
                "read_name": read.query_name,
                "chrom": chrom,
                "start": start,
                "end": end,
                "orientation": orientation,
                "insert_size": abs(read.template_length)
            })

    bam.close()

    if len(records) == 0:

        df_discordant_reads_BED = pd.DataFrame(columns=['chrom', 'start', 'end'])

    else:

        df = pd.DataFrame(records)
        
        df_discordant_reads_BED = df[['chrom', 'start', 'end']]
        
        # convert to 0-indexed half open intervals
        df_discordant_reads_BED['start'] -= 1
        df_discordant_reads_BED[['start', 'end']] = df_discordant_reads_BED[['start', 'end']].astype(int)
    
    # save BED file
    df_discordant_reads_BED.to_csv(discordant_alns_BED_file, sep='\t', header=None, index=False)
    
    # run bedtools genomecov
    ref_genome_chrom_lengths_file = f"{ref_genome_file.replace('.fna', '')}.chrom_lengths.txt"
    
    # then run the following to compute the number of discordant reads
    if not os.path.isfile(ref_genome_chrom_lengths_file):
        subprocess.run(f"cut -f1,2 {ref_genome_file}.fai > {ref_genome_chrom_lengths_file}", shell=True)

    subprocess.run(f"bedtools genomecov -i {discordant_alns_BED_file} -g {ref_genome_chrom_lengths_file} -d | gzip -c > {discordant_alns_by_pos_file}", shell=True)
    


############################ GOAL: TO GET 3 ALIGNMENT STATISTICS FROM THE BAM FILE ############################


############################ STEP 1: ROLLING AVERAGE OF COVERAGE FROM LEFT AND RIGHT SIDES ############################


#depth_file = os.path.join(os.path.dirname(bam_file), f"{sample}.depth.tsv.gz")
depth_file = f"/n/scratch/users/s/sm624/FP_characteristics/{simulator}_{genome}/bam_depths/{sample}.depth.tsv.gz"
df_depth = pd.read_csv(depth_file, compression='gzip', sep='\t', header=None, names=['CHROM', 'POS', 'COV'])

# compute rolling average of coverage. Smaller than 100 is too small, not enough smoothing
window_size = 100
df_depth['COV_LEFT_ROLLING_AVG'] = df_depth['COV'].rolling(window=window_size, min_periods=1, closed='right').mean()

# for this one, compute the rolling average as you would from the left, but reverse the values beforehand. This is the easiest way. Then reverse them again
df_depth['COV_RIGHT_ROLLING_AVG'] = df_depth['COV'][::-1].rolling(window=window_size, min_periods=1, closed='right').mean()[::-1]

# take the maximum at each site
df_depth['COV_MAX_ROLLING_AVG'] = np.max([df_depth['COV_LEFT_ROLLING_AVG'], df_depth['COV_RIGHT_ROLLING_AVG']], axis=0)

    
############################ STEP 2: NUMBER OF SOFT-CLIPPED BASES ADJACENT TO OR AT EACH VARIANT SITE ############################


# generate the file of soft-clipped bases (analogous to the depth file, with number of soft-clipped bases at each site
if not os.path.isfile(soft_clips_by_pos_file):
    save_table_of_soft_clips(sample, bam_file, genome_file)

df_soft_clips = pd.read_csv(soft_clips_by_pos_file, compression='gzip', sep='\t', header=None, names=['CHROM', 'POS', 'CLIPPED_BASES'])


############################ STEP 3: NUMBER OF DISCORDANTLY PAIRED READS AT EACH SITE ############################


# generate the file of soft-clipped bases (analogous to the depth file, with number of soft-clipped bases at each site
if not os.path.isfile(discordant_alns_by_pos_file):
    save_table_of_discordant_reads(sample, bam_file, genome_file)

df_discordant_alns = pd.read_csv(discordant_alns_by_pos_file, compression='gzip', sep='\t', header=None, names=['CHROM', 'POS', 'DISCORDANT_READS'])


############################ STEP 4: PUT IT ALL TOGETHER ############################

## merge the coverage, soft clipped bases and discordant reads information
df_depth.set_index("POS", inplace=True)
df_soft_clips.set_index("POS", inplace=True)
df_discordant_alns.set_index("POS", inplace=True)

df_genome_metrics = df_depth[['COV', 'COV_MAX_ROLLING_AVG']].merge(df_soft_clips[['CLIPPED_BASES']],  how='left', left_index=True, right_index=True).merge(df_discordant_alns[['DISCORDANT_READS']], how='left', left_index=True, right_index=True)

# compute the ratio
df_genome_metrics['COV_RATIO'] = df_genome_metrics['COV'] / df_genome_metrics['COV_MAX_ROLLING_AVG']

df_genome_metrics['COV'] = df_genome_metrics['COV'].replace(0, np.nan) #to prevent ZeroDivisionError

# normalize to coverage so that differences in coverage between samples doesn't swamp the signal
df_genome_metrics['DISCORDANT_READS_RATIO'] = df_genome_metrics['DISCORDANT_READS'] / df_genome_metrics['COV']
df_genome_metrics['CLIPPED_BASES_RATIO'] = df_genome_metrics['CLIPPED_BASES'] / df_genome_metrics['COV']

# write final file
df_genome_metrics.to_csv(f"{out_dir}/{sample}.genome_metrics.csv", index=True)