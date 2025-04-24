# %%
import pysam
import os
import pandas as pd
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# %%
path="path to results/"
result_df= pd.read_csv(path+"IS_frequencies.tab",sep="\t")
df_IS = pd.read_csv(path+"all_IS.tab",sep="\t")

# %%
def is_partially_mapped(cigar):
    """
    Determines if a sequencing read is partially mapped based on its CIGAR string.

    A read is considered partially mapped (soft-clipped) if either the start or the end of the read
    has more than 20 bases that are not aligned to the reference genome. This typically indicates
    that only a portion of the read aligns to the reference, suggesting potential structural variations
    or insertions/deletions at the alignment location.

    Parameters:
    - cigar (list of tuples): The CIGAR string of the read, parsed into a list of tuples where each tuple
      contains an operation code (int) and the length (int) of that operation.

    Returns:
    - bool: True if the read is partially mapped (soft-clipped) by more than 20 bases at either end,
            False otherwise.
    """
    # CIGAR operation codes: 4 represents 'Soft-clipped' in pysam and SAM/BAM file formats.
    if (cigar[0][0] == 4 and cigar[0][1] > 20) or (cigar[-1][0] == 4 and cigar[-1][1] > 20):
        return True
    return False
    
def get_heterozygosity(path,sample:str,variant, pos:str, contig, start, end, mask=False, min_ratio=0.3, TSD=10, mapq=5):
    """
    Analyzes sequencing data to determine the zygosity of a specified variant within a given genomic region.

    This function extracts reads aligning to a specified region, then analyzes these reads to determine
    if the variant is heterozygous, homozygous for the insertion, or homozygous with no insertion, based
    on the presence and distribution of reads crossing or clipped near the variant position.

    Parameters:
    - path (str): Base directory path for SAM/BAM file storage (not directly used in the snippet).
    - sample (str): Identifier for the sample being analyzed.
    - variant (str): Identifier for the variant.
    - pos (int or str): Position of the variant. If a range (str with '-'), it indicates a region.
    - contig (str): The name of the chromosome or contig where the variant is located.
    - start (int): Start position of the genomic region of interest.
    - end (int): End position of the genomic region of interest.
    - mask (bool): Whether specific regions in the assembly are masked to viral sequences.
    - min_ratio (float): Minimum ratio of crossing reads to consider the variant heterozygous.
    - TSD (int): Size threshold for target site duplication (TSD) to consider in analysis.
    - mapq (int): Minimum mapping quality of reads to include in the analysis.

    Returns:
    - int: An integer code representing the zygosity determination:
        - 2: Homozygous for the insertion.
        - 1: Heterozygous.
        - -1: No data or unable to determine.
        - -2: Homozygous with no insertion, or very likely a somatic event.
    """
    try:
        if pos!=None and pos!='nan' and pos!='' and pos!='NAN' and pos!='NaN':
            os.system("samtools view -h path_to_alignments/"+str(sample)+"_sorted.bam "+contig+":"+str(start)+"-"+str(end)+" > "
            +path+sample+"_"+variant+"_"+contig+":"+str(start)+"-"+str(end)+"_"+sample+".sam")
            samfile = pysam.AlignmentFile(path+sample+"_"+variant+"_"+contig+":"+str(start)+"-"+str(end)+"_"+sample+".sam", "r")
            data = []

            # Iterate through each alignment
            for read in samfile.fetch():
                if not read.is_unmapped:
                    # Extracting data from the SAM file and mapping to desired columns
                    qseq = read.query_name
                    qlen = read.query_length
                    qstart = read.query_alignment_start
                    qend = read.query_alignment_end
                    tstrand = '-' if read.is_reverse else '+'
                    tseq = samfile.get_reference_name(read.rname)
                    tlen = read.reference_length
                    tstart = read.reference_start
                    tend = read.reference_end
                    bases = read.query_alignment_sequence
                    mapq = read.mapping_quality
                    NM = read.get_tag('NM') if read.has_tag('NM') else None
                    AS = read.get_tag('AS') if read.has_tag('AS') else None
                    soft= is_partially_mapped(read.cigar)
                    data.append([qseq, qlen, qstart, qend, tstrand, tseq, tlen, tstart, tend, bases, mapq, NM, AS,soft])

            df = pd.DataFrame(data, columns=['qseq', 'qlen', 'qstart', 'qend','tstrand', 'tseq', 'tlen', \
                                                'tstart', 'tend', 'bases', 'mapq', 'NM', 'AS','soft'])

            df=df[df['mapq']>=mapq].reset_index(drop=True)
            samfile.close()
            os.system("rm "+path+sample+"_"+variant+"_"+contig+":"+str(start)+"-"+str(end)+"_"+sample+".sam")
            if len(df)==0:
                return -1

            # The function includes logic to handle different scenarios based on the 'pos' parameter,
            # whether it specifies a single position or a range, and whether to consider reads as supporting
            # heterozygosity or homozygosity based on their alignment characteristics and the specified thresholds.
            if '-' in pos and mask: #That means there is a long gap between the breakpoints, due to the insertion also being in the assembly, we check if any reads are clipped in the breakpoints
                pos1=int(pos.split("-")[0]);pos2=int(pos.split("-")[1])
                clipped = len(df[((abs(df['tstart']-pos1)<=TSD) | (abs(df['tend']-pos1)<=TSD)) & df['soft']])
                crossing = len(df[((df['tstart']<(pos1-TSD)) & (df['tend']>(pos1+TSD))) & ~df['soft']])
                total = crossing + clipped
                ratio = crossing/total
                if ratio>(1-min_ratio):
                    return 2 #homozygous insertion
                elif ratio>min_ratio:
                    return 1 #heterozygous
                else: # try checking pos2
                    clipped = len(df[((abs(df['tstart']-pos2)<=TSD) | (abs(df['tend']-pos2)<=TSD)) & df['soft']])
                    crossing = len(df[((df['tstart']<(pos2-TSD)) & (df['tend']>(pos2+TSD))) & ~df['soft']])
                    total = crossing + clipped
                    ratio = crossing/total
                    if ratio>(1-min_ratio):
                        return 2 #homozygous insertion
                    elif ratio>min_ratio:
                        return 1 #heterozygous
                    else:
                        return -2 #homozygous no insertion --> not able to confirm the insertion (very likely a somatic, code -2 to count later on somatic events

            else: #that means this insertion is not in the assembly, and we check if any reads cross the insertion site
                if '-' in pos:
                    pos = int(pos.split("-")[0])
                else:
                    pos=int(pos)
                clipped = len(df[((abs(df['tstart']-pos)<=TSD) | (abs(df['tend']-pos)<=TSD)) & df['soft']])
                crossing = len(df[((df['tstart']<(pos-TSD)) & (df['tend']>(pos+TSD))) & ~df['soft']])
                total = crossing + clipped
                if total==0:
                    return -1
                ratio = crossing/total
                if ratio<min_ratio:
                    return 2 #homozygous insertion
                elif ratio>(1-min_ratio):
                    return -2 #homozygous no insertion --> not able to confirm the insertion (very likely a somatic, code -2 to count later on somatic events)
                else:
                    return 1 #heterozygous          
        return -1  # Default return value if conditions for zygosity determination are not met.
    except Exception as e:
        # Handle or log the exception
        # Consider returning a simple error message or a default value
        print(f"Error: {e}")
        return None  # Or some error indicator
path="path to processed data/"
for i,group in enumerate(result_df['group']):
    args=[]
    temp_samples = df_IS[df_IS['group']==group]['Sample'].unique()
    for sample in temp_samples:
        temp = df_IS[(df_IS['group']==group) & (df_IS['Sample']==sample)].reset_index(drop=True)
        args.append((path,str(sample),str(result_df.loc[i,"integration_variant"]),str(temp.loc[0,'POS']),str(temp.loc[0,'tseq']),int(temp.loc[0,'tstart']),int(temp.loc[0,'tend']),temp.loc[0,'masked']))
    with Pool(15) as p:
        zygosity = p.starmap_async(get_heterozygosity,args).get()
    for j,sample in enumerate(temp_samples):
        result_df.loc[i,sample]=zygosity[j]
samples = df_IS['Sample'].unique() #ignoring the duplicates at this point
result_df['#Missing']= result_df[samples].apply(lambda row: (row == -1).sum(), axis=1)
result_df['#Somatic']= result_df[samples].apply(lambda row: (row == -2).sum(), axis=1)
result_df['AAF']= result_df[samples].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
result_df.to_csv(path+"IS_frequenciesv2.tab",sep="\t",index=False)


