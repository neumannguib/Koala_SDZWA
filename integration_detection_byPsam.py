import pysam
import glob
import pandas as pd
import os
from collections import Counter
from itertools import combinations
from scipy.stats import mode

threads = 6
min_reads = 20 ## amount of reads to provide evidence for KoRV insertion threshold
global pos_threshold_others
global pos_threshold_korv
pos_threshold_korv=500 ## threshold for distance of reads aligned in genome to identify KoRV insertion site (radius like; because of paired-short-reads)
pos_threshold_others= 9000
mapq_filter = 30


#general function to parse sam file
def sam_parse(sam,mapq_filterthreshold = 0):
    """ Parse sam file to a pandas dataframe, filtering based on mapping quality

    :param sam: str, path to sam file
    :param mapq_filterthreshold: int, filter threshold of mapping quality
    :return: pd.DataFrame, filtered alignments 
    """
    
    samfile = pysam.AlignmentFile(sam, "r")
    data = []

    # Iterate through each alignment
    for read in samfile.fetch():
        if not read.is_unmapped:
            # Extracting data from the SAM file and mapping to desired columns
            qseq = read.query_name
            qlen = read.query_length
            qstart = read.query_alignment_start
            qend = read.query_alignment_end
            strand = '-' if '_R' in qseq else '+'
            read1 = 1 if '_Read1' in qseq else 2
            tstrand = '-' if read.is_reverse else '+'
            soft = True if "_softclipped" in qseq else False
            side = 'right' if 'right' in qseq else 'left' if 'left' in qseq else None
            tseq = samfile.get_reference_name(read.rname)
            tlen = read.reference_length
            tstart = read.reference_start
            tend = read.reference_end
            bases = read.query_alignment_sequence
            mapq = read.mapping_quality
            variant = qseq.split("_")[-2]
            gene = qseq.split("-")[-1] if 'soft' in qseq else ""
            region = qseq.split("_")[-1].split("-")[0]+"-"+qseq.split("_")[-1].split("-")[1]
            NM = read.get_tag('NM') if read.has_tag('NM') else None
            AS = read.get_tag('AS') if read.has_tag('AS') else None

            data.append([qseq, qlen, qstart, qend, strand, read1, tstrand, soft, side, tseq, tlen, tstart, tend, bases, mapq, variant, region, gene, NM, AS])

    df = pd.DataFrame(data, columns=['qseq', 'qlen', 'qstart', 'qend', 'strand', 'read_in_pair','tstrand', 'softclipped','side','tseq', 'tlen', \
                                 'tstart', 'tend', 'bases', 'mapq',"integration_variant","Viral_region","Viralgene", 'NM', 'AS'])

    samfile.close()
    df = df[df.mapq >= mapq_filterthreshold].reset_index() ## Filter on mapq threshold
    return df


def find_multiple_modes(numbers, min_difference=10, min_occurrences=4):
    # Count occurrences of each number
    counts = Counter(numbers)
    
    # Find all values with occurrences >= min_occurrences
    valid_values = [num for num, count in counts.items() if count >= min_occurrences]
    
    # Sort valid_values based on their counts, from most to least common
    valid_values_sorted = sorted(valid_values, key=lambda x: counts[x], reverse=True)
    
    # Check combinations of valid values for a difference greater than min_difference
    valid_modes = []
    for val1, val2 in combinations(valid_values_sorted, 2):
        if abs(val1 - val2) > min_difference and val1 not in valid_modes and val2 not in valid_modes:
            valid_modes.extend([val1, val2])
    
    # Remove duplicates and sort by abundance
    valid_modes = sorted(list(set(valid_modes)), key=lambda x: counts[x], reverse=True)
    
    return valid_modes

def identify_integrationsites(sam_df: pd.DataFrame, min_reads: int):
    """finds all insertions with at least <min_reads> providing evidence for given insertion

    :param sam_df: pd.DataFrame, .sam file dataframe, by sam_parse() 
    :param min_reads: int, amount of reads to provide evidence for KoRV insertion threshold
    :param pos_threshold: int, distance threshold to cluster together reads
    :return clustered_reads: pd.DataFrame of identified integration sites (positions just most common start and end positions of alignments)
    """
    def assign_cluster(group):
        # Get the pos_threshold for the current integration_variant
        variant = group['integration_variant'].iloc[0]  # Assuming all values are the same within the group
        #pos_threshold = pos_threshold_korv if variant=='KoRV' else pos_threshold_korv if variant=='PhaCinBeta' else pos_threshold_others
        pos_threshold = pos_threshold_korv if 'PhaCinBeta' not in variant else pos_threshold_others
        
        # Calculate diff, compare with pos_threshold, and cumsum for cluster assignment
        cluster_ids = group['tstart'].diff().gt(pos_threshold).cumsum()
        return cluster_ids

    # Apply the custom function to each group and assign the result to a new column
    sam_df['reads_cluster'] = sam_df.sort_values(["integration_variant", "tseq", "tstart"])\
                                .groupby(["tseq", "integration_variant"], group_keys=False)\
                                .apply(assign_cluster)

    clustered_reads = []
    ## Loop through all clusters of reads with k = (contigname, clusterID) and v = [indices of rows in this cluster, ]
    for k, v in sam_df.groupby(['tseq', 'integration_variant','reads_cluster']).groups.items():
        reads = len(v)
        if reads >= min_reads:
            mode_1 =sam_df.loc[sam_df.index.isin(v) & (((sam_df['side']=='right') & (sam_df['tstrand']=='+')) | 
           ((sam_df['side']=='left') & (sam_df['tstrand']=='-')) ), 'tstart']
            mode_2 =sam_df.loc[sam_df.index.isin(v) & (((sam_df['side']=='left') &  (sam_df['tstrand']=='+')) | 
            ((sam_df['side']=='right') & (sam_df['tstrand']=='-'))), 'tend']
            n=len(mode_1)+len(mode_2)
            if n>1:
                # Find cases of multiple modes 500bp away 
                valid_modes = find_multiple_modes(list(mode_1)+list(mode_2))
                if len(valid_modes)==0:
                    modes=mode(list(mode_1)+list(mode_2))[0][0] #get mode   (novel insertion)
                else:
                    multiple=sorted([valid_modes[0],valid_modes[1]]) #get breakpoints (case when reference genome already contains the viral insertion)
                    modes=str(multiple[0])+"-"+str(multiple[1])
            else:
                modes = None
                n=0
            clustered_reads.append([k[0], 
                                    modes,
                                    n,
                                    sam_df.iloc[v]['tstart'].min(), 
                                    sam_df.iloc[v]['tend'].max(), 
                                    k[1], reads,str(k[2])+"_"+k[0]+"_"+k[1]])
    clustered_reads = pd.DataFrame(clustered_reads, columns=['tseq','POS', '#clipped_reads','tstart', 'tend', 'integration_variant', 
    'reads_evidence','reads_cluster'])
    return clustered_reads

def map_2_koala(file): 
    path = "path to processed data"
    sample = str(file.split('/')[-1].split('_IS')[0])
    ## process alignment to Koala 
    sam_df = sam_parse(path+sample+"_koala_mapped.sam", mapq_filter)

    ## Identify KoRV integration sites
    #sam_df['integration_variant']=sam_df['integration_variant'].replace({"AB721500.1":"KoRV","AF151794.2":"KoRV","KC779547.1":"KoRV"})
    sam_df['integration_variant']=sam_df['integration_variant'].replace({"gi|425029687|dbj|AB721500.1|":"KoRV-A",
    "gi|516464655|gb|KC779547.1|":"KoRV-B"})
    sam_df['integration_variant']=sam_df['integration_variant'].apply(lambda x: "PhaCinBetalike" if "Betalike" in x else "PhaCinBeta" if "Beta" in x else x)
    integration_sites = identify_integrationsites(sam_df, min_reads)
    is_tsv = path+sample+"_integration_sites.tsv"
    integration_sites.to_csv(is_tsv, sep='\t', index=False)

path = "path to processed data"
names = sorted(glob.glob(path+"*.fastq"))
from multiprocessing import Pool
with Pool(5) as p:
    p.map(map_2_koala,names)