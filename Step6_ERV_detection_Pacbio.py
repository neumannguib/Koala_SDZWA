# %% [markdown]
# **Pipeline KoRV integration site discovery from PacBio reads**
# Modified from Max Coenen, by Guilherme Neumann

# %%
import glob
import pandas as pd
import gzip
from Bio import SeqIO
from functools import reduce
import numpy as np
from subprocess import check_output
from subprocess import DEVNULL
import os
from multiprocessing import Pool
import warnings
warnings.filterwarnings('ignore')

# %%
### KoRV variants & phaCin-Beta
korv_fa = "Retroviruses/KoRV_phaCinBeta.fa"
korv_db_path = "Retroviruses/blastdb_KoRV_phaCinBeta"
 # makeblastdb -in Retroviruses/KoRV_phaCinBeta.fa -dbtype nucl 
 # -out Retroviruses/blastdb_KoRV_phaCinBeta

### koala reference genome
koala_ref = "GCA_030178435.1_masked_KoRV_PhaCinBeta/GCA_030178435.1_ASM3017843v1_genomic.fna.masked"

primers = "primers.fa"
threads = 25
min_reads = 10 ## amount of reads to provide evidence for KoRV insertion threshold per start and end positions
len_match = 50 #minimum blast match in bp
identity = 90 #minimum identify of blast results

### OUTPUT DIRECTORY:
out_dir = "data/processed/"

## Reads
fastqs = glob.glob("data/raw/*.fastq.gz") 

# %% [markdown]
# ### Pipeline that starts from reads, and identifies KoRV integration sites ###

# %%
def read_slicer(grp_blast: pd.core.groupby.generic.DataFrameGroupBy, record_dict: dict):
    """Takes largest unaligned (to KoRV) part of read, and creates fastq format of that softclipped part

    :param grp_blast: pd.core.groupby.generic.DataFrameGroupBy, applied on pd.DataFrame, of blastn format 6 in as df
    :param record_dict: SeqIO.RecordDict
    :return fastq: string, 
    """
    ## Identify best alignment to KoRV, and set `grp_blast` to respective target sequence alignment(s)
    korv = grp_blast.loc[(grp_blast.saccver != "PhER") & (grp_blast.pident>=identity)]
    if korv.empty:
        return ''
    best_al = korv.loc[korv.sort_values('saccver')['bitscore'].idxmax()]
    ## check if read also aligns to PhER sequence
    pher = grp_blast.loc[(grp_blast.saccver == "PhER") & (grp_blast.pident>=identity)]
    if not pher.empty:
        best_al.saccver += "PhER"

    ## indentify parts of read that can be aligned to KoRV ref. sequence, and part of read to map to genome
    ranges = tuple(zip(grp_blast.qstart, grp_blast.qend))
    aligned = set()
    [aligned.update(tuple(range(al[0], al[1] + 1))) for al in ranges]

    aligned = set(aligned) 
    unaligned = set(range(1, best_al.qlen + 1)) - aligned
    if unaligned == set(): ## When the entire read aligns to KoRV & PhER -> no flanking sequence
        return ''

    ## Identify longest part of read that is unaligned to KoRV
    u_ranges, un_i = [], -1
    for v in sorted(unaligned):
        if v - 1 not in unaligned:
            un_i += 1
            u_ranges.append( [v, 0] )
        else:
            u_ranges[un_i][1] = v
    longest_range = sorted([((therange[1] - therange[0]), therange) for therange in u_ranges], reverse = True)[0]    
    
    ## Construct fastq of part of read unaligned to KoRV
    if longest_range[0] >= 20:
        softclip_seq = record_dict[best_al.qaccver].seq[longest_range[1][0] -1: longest_range[1][1]]
        softclip_q = ''.join(list(map(lambda x: chr(x+33), 
                                      record_dict[best_al.qaccver].letter_annotations["phred_quality"][longest_range[1][0]-1: longest_range[1][1]])))
        #side = 'left' if (best_al.qlen-longest_range[1][1])<len_match else 'rigth'
        ## Header: @Read_id/lengthofOGread/posinOGread/KoRVvariant     
        fastq_rec = f"@{best_al.qaccver}/len{best_al.qlen}/slicepos{longest_range[1][0]}-{longest_range[1][1]}/{best_al.saccver}\n{softclip_seq}\n+\n{softclip_q}\n"
        return fastq_rec
    else:
        return ''


# %%
def paf_parse(paf: str):
    """Parses .paf file and 

    :param paf: str, path to .paf file (minimap2)
    :param min_reads: int, amount of reads to provide evidence for KoRV insertion threshold
    :return insertions: pd.DataFrame of (KoRV) insertions in genome supported by <min_reads> amount of reads 
    """
    paf_df = pd.read_csv(paf, sep='\t', header=None, \
                         names = ['qseq', 'qlen', 'qstart', 'qend', 'strand', 'tseq', 'tlen', \
                                'tstart', 'tend', 'match', 'bases', 'mapq', 'NM', 'ms', 'AS', \
                                'nn', 'tp', 'cm', 's1', 's2', 'de', 'rl', 'cg'], \
                         usecols = ['qseq', 'qlen', 'qstart', 'qend', 'strand', 'tseq', 'tlen', \
                                  'tstart', 'tend', 'match', 'bases', 'mapq', "AS"])
    paf_df["AS"] = paf_df["AS"].str.replace("AS:i:", "").astype(int) ## Alignment score tag
    paf_df["integration_variant"] = paf_df['qseq'].str.rpartition('/')[2]
    ## Additionally, the r1 tag might be interesting: Length of query regions harboring repetitive seeds    
    return paf_df

# %%
def identify_integrationsites(paf_df: pd.DataFrame, min_reads: int):
    """finds all insertions with at least <min_reads> providing evidence for given insertion

    :param paf_df: pd.DataFrame, .paf file dataframe, by paf_parse() (minimap2)
    :param min_reads: int, amount of reads to provide evidence for KoRV insertion threshold
    """
    from scipy.stats import mode
    def grouped_mode(ser):
        mod, count = mode(ser)
        if isinstance(mod, np.ndarray):
            return mod[0] if len(mod) > 0 else None
        else:
            return mod  # In case mod is a scalar

    def grouped_mode_count(ser):
        mod, count = mode(ser)
        if isinstance(count, np.ndarray):
            return count[0] if len(count) > 0 else None
        else:
            return count  # In case mod is a scalar
    ## Select insertions with at least min_reads amount of reads at that position
    # Here, we will create bins with a certain range width (10 in this case)
    range_width = 10 # allow reads to be 10bp of distance in the start or end position
    #filter based on mapquality
    paf_df = paf_df[(paf_df["mapq"]>=30) & (paf_df["bases"]>=len_match) & ((paf_df["match"]/paf_df["bases"])>=0.9)]
    paf_df['bin'] = (paf_df['tstart'] // range_width) * range_width
    grouped = paf_df.groupby(['integration_variant', 'tseq', 'bin']).filter(lambda x: len(x) >= min_reads)
    tstart = grouped.groupby(['integration_variant', 'tseq', 'bin']).agg(tstart=('tstart', grouped_mode),
    block_start = ('tstart', 'min') ,tstart_reads_evidence=('tstart', 'size'), coverage_tstart=('tstart', grouped_mode_count)).reset_index()

    paf_df['bin'] = (paf_df['tend'] // range_width) * range_width
    grouped = paf_df.groupby(['integration_variant', 'tseq', 'bin']).filter(lambda x: len(x) >= min_reads)
    tend = grouped.groupby(['integration_variant', 'tseq', 'bin']).agg(tend=('tend', grouped_mode), block_end = ('tend', 'max'),
    tend_reads_evidence=('tend', 'size'), coverage_tend=('tend', grouped_mode_count)).reset_index()

    ## Find KoRV integration sites that have evidence from reads on the 5'- and 3'-side 
    #insertions = tstart.merge(tend, on=["tseq", "integration_variant"], copy=False).loc[lambda x: abs(x['tstart'] - x['tend']) <= 10]
    insertions = tstart.merge(tend, on=["tseq", "integration_variant"], copy=False).loc[lambda x: abs(x['tstart'] - x['tend']) <= 1411]#1411 longest recKoRV masked
    ## Format layout of dataframe
    insertions.rename(columns = {"tseq_x":'tseq'}, inplace = True)
    #insertions.rename(columns = {"0_x":'tstart_reads_evidence'}, inplace = True)
    #insertions.rename(columns = {"0_y":'tend_reads_evidence'}, inplace = True)
    ## Filter such that integration sites that have alternative tend for a given tstart are removed based on highest number of evidence reads
    insertions = insertions.groupby(['integration_variant', 'tseq', 'tstart'], as_index=False).apply(lambda x: x.loc[x.tend_reads_evidence.idxmax()])
    insertions = insertions.groupby(['integration_variant', 'tseq', 'tend'], as_index=False).apply(lambda x: x.loc[x.tstart_reads_evidence.idxmax()])
    insertions['coverage'] = insertions['tstart_reads_evidence'] +  insertions['tend_reads_evidence']
    
    # Sort the DataFrame by tseq, tstart, and descending coverage
    sorted_df = insertions.sort_values(by=['tseq', 'tstart', 'coverage'], ascending=[True, True, False])
    filtered_insertions = []
    # remove overlapping insertions
    for tseq, group in sorted_df.groupby('tseq'):
        # Initialize a variable to keep track of the end of the last added interval
        last_end = -1
        for _, row in group.iterrows():
            # If the current tstart is greater than the end of the last added interval, add the row
            if row['tstart'] > last_end:
                filtered_insertions.append(row)
                # Update last_end to the end of the current interval
                last_end = row['tend']
    return pd.DataFrame(filtered_insertions)

# %%
## Main, perform pipeline for all samples
def main(reads):
    sample_no = reads.split('/')[-1].split('.')[0]

    print(f"Processing (trimming & BLAST to KoRV): {reads}")
    trimmed_reads = out_dir + 'trimmed_'+reads.split('/')[-1].split('.')[0]

    ## trim reads - just removing optical duplications
    os.system("fastp -i "+reads+" --stdout --dedup --thread "+str(threads)+" \
    --json "+trimmed_reads+".json --html "+trimmed_reads+".html > "+trimmed_reads+".fastq")
    os.system("cat "+trimmed_reads+".fastq | sed -n '1~4s/^@/>/p;2~4p' >"+trimmed_reads+".fasta")

    ## BLAST fastq reads against KoRV & PHER reference
    
    blastout = f"{out_dir}KoRV_Blast/{sample_no}_korv_blast.txt"

    if not os.path.isdir(f"{out_dir}KoRV_Blast"): 
        check_output(f"mkdir {out_dir+'KoRV_Blast'}", shell=True, stderr=DEVNULL)
    blast_cmd= f'blastn -query {trimmed_reads}.fasta -perc_identity 90 -mt_mode 1 -subject_besthit \
                    -out {blastout} -db {korv_db_path} \
                    -outfmt "6 qaccver saccver pident length mismatch gapopen qlen qstart qend sstart send sstrand evalue bitscore" \
                    -num_threads {threads}'
    check_output(blast_cmd, shell=True, stderr=DEVNULL)
    ## reads blast results into pd.Dataframe
    blast = pd.read_csv(blastout, sep='\t', names=["qaccver","saccver","pident","length","mismatch","gapopen",
    "qlen","qstart","qend","sstart","send","sstrand","evalue","bitscore"])
    ## dictionary of entire fastq file (RAM on the server easily allows this)
    record_dict = SeqIO.to_dict(SeqIO.parse(trimmed_reads+".fastq", "fastq"))
    ## create fastq of softclipped reads
    softclip_fastq = blast.groupby('qaccver').apply(lambda x: read_slicer(x, record_dict))  

    softclip_fastq_path = f"{out_dir}softclipped_fastqs/{sample_no}_softclip.fastq"
    if not os.path.isdir(f"{out_dir}softclipped_fastqs"): 
        check_output(f"mkdir {out_dir+'softclipped_fastqs'}", shell=True, stderr=DEVNULL)
    with open(softclip_fastq_path, 'w') as fq_out:
        fq_out.writelines(softclip_fastq)

    
    ## Align softclipped parts to Koala Genome
    paf_out = f"{out_dir}genome_paf/{sample_no}_koalaKoRV.paf"
    if not os.path.isdir(f"{out_dir}genome_paf"): 
        check_output(f"mkdir {out_dir+'genome_paf'}", shell=True, stderr=DEVNULL)
    mm_cmd = f"minimap2 -x map-pb -Y -c --secondary=no  -t {threads} {koala_ref} {softclip_fastq_path} > {paf_out}"
    check_output(mm_cmd, shell=True, stderr=DEVNULL)
    
    ## process alignment to Koala (paf)
    paf_df = paf_parse(paf_out)

    ## Identify KoRV integration sites
    integration_sites = identify_integrationsites(paf_df, min_reads)
    is_tsv = f"{out_dir}{sample_no}_integration_sites.tsv"
    integration_sites.to_csv(is_tsv, sep='\t', index=False)
with Pool(10) as p:
    p.map(main,fastqs)

# %% [markdown]
# ### Comparing results between samples ###

# %%
def check_overlap_results(df1: pd.DataFrame, sample1: str, df2: pd.DataFrame, sample2: str):
    """Compares results between two samples, and creates pd.DataFrame where all identical integration sites between samples are stored

    :param df1: pd.DataFrame, pf.Dataframe of integration sites by identify_integrationsites()
    :param sample1: string, name of sample 1
    :param df2: pd.DataFrame, pf.Dataframe of integration sites by identify_integrationsites()
    :param sample2: string, name of sample 2
    """
    overlap = df1.merge(df2, how='inner', on=['integration_variant', 'tseq', 'tstart', 'tend'])
    overlap.rename(columns = {"tstart_reads_evidence_x":f'tstart_reads_evidence_{sample1}', 
                              "tstart_reads_evidence_y":f'tend_reads_evidence_{sample2}',
                              "tend_reads_evidence_x":f'tend_reads_evidence_{sample1}', 
                              "tend_reads_evidence_y":f'tend_reads_evidence_{sample2}'}, inplace = True)
    return overlap

# %%
def compare_results(all_results: dict):
    """ Compares results between all samples, and creates pd.DataFrame where values represent pairwise amount of identical integration sites between samples

    :param all_results: dict, {'sample_name': pf.Dataframe of integration sites by identify_integrationsites()}
    :return overlap_df: pd.DataFrame where columns and indices are samples, and values are pairwise identical integration sites
    """
    overlap_df = pd.DataFrame(columns=list(all_results.keys()), index=list(all_results.keys()))
    for sample, r_df in all_results.items():
        for sample2, r_df2 in all_results.items():
            ## See what the overlap is between two samples 
            overlap = check_overlap_results(r_df, sample, r_df2, sample2)
            overlap_df[sample][sample2] = len(overlap)
            overlap_df[sample2][sample] = len(overlap)
    return overlap_df

# %%
## Reads all integration_sites.tsv files
results = sorted(glob.glob(f"{out_dir}*_integration_sites.tsv"))

all_results = {};df_IS = pd.DataFrame()
for result in results:
    if "Milla" not in result and "Bunji" not in result:
        sample_no = result.split('/')[-1].split('_')[1]
        r_df = pd.read_csv(result, sep='\t')
        all_results[sample_no] = r_df
        r_df['Sample'] = sample_no
        df_IS = pd.concat([df_IS,r_df])
all_integration_sites = compare_results(all_results)
all_is_tsv = f"shared_integration_sites.tsv"
all_integration_sites.to_csv(all_is_tsv, sep='\t')
all_integration_sites

# %%
len(df_IS['Sample'].unique()) #number of koala samples

# %%
df_IS=df_IS.reset_index(drop=True)
df_IS.groupby('integration_variant').count()['tstart'] #total number of IS per type for all samples, KC779547.1 KoRVB

# %%
df_IS=df_IS.reset_index(drop=True)
df_IS.groupby('integration_variant').count()['tstart'] #total number of IS per type for all samples, KC779547.1 KoRVB

# %%
df_IS['windown_len']=abs(df_IS['tend']-df_IS['tstart'])
df_IS['block_len']=abs(df_IS['block_end']-df_IS['block_start'])
df_IS.describe()

# %%
# Sort the DataFrame by POSITIONS + TSD
df_IS['integration_variant']=df_IS['integration_variant'].replace({'gi|425029687|dbj|AB721500.1|':'KoRV-A','gi|516464655|gb|KC779547.1|':'KoRV-B'})
TSD=10
df_IS["strand"]=[True if df_IS.loc[i,'tstart']<df_IS.loc[i,'tend'] else False for i in df_IS.index]
df_IS['start']=[df_IS.loc[i,'tstart'] if df_IS.loc[i,'tstart']<df_IS.loc[i,'tend'] else df_IS.loc[i,'tend'] for i in df_IS.index]
df_IS['end']=[df_IS.loc[i,'tend'] if df_IS.loc[i,'tstart']<df_IS.loc[i,'tend'] else df_IS.loc[i,'tstart'] for i in df_IS.index]
df_IS = df_IS.sort_values(by=['integration_variant','tseq','start','end']).reset_index(drop=True)

# Initialize a column to keep track of merged groups
df_IS['group'] = range(1, len(df_IS) + 1)  # Each row starts in its own group

# Function to find overlapping groups and update the 'group' column based on tstart and tend of clusters 
def update_groups(row):
    current_group = row['group']
    temp = df_IS[(df_IS['tseq']==row['tseq']) & (df_IS['integration_variant']==row['integration_variant']) & (df_IS['end'] >= row['start']) & 
    (df_IS['start'] <= row['end']) &
        (df_IS['group'] != current_group)]
    overlapping_groups = temp['group'].tolist()
    if overlapping_groups:
        min_group = min([current_group] + overlapping_groups)
        df_IS.loc[df_IS['group'].isin([current_group] + overlapping_groups), 'group'] = min_group

# Iterate through rows and update the 'group' column
df_IS.apply(update_groups, axis=1)
free_groups=[group for group in range(1, len(df_IS) + 1) if group not in df_IS["group"].unique()]
#in case one sample has more than one IS overlapping to others, I keep only the one the highest coverage, the most likely ERV
for i in df_IS[df_IS.sort_values("tstart_reads_evidence",ascending=False).duplicated(["Sample","group"])].index:
    df_IS.loc[i,"group"]=free_groups.pop()
    
def get_mode(series):
    modes = series.mode()
    if not modes.empty:
        return modes.iloc[0]
    else:
        return None
def get_mode_count(series):
    modes = series.mode()
    if not modes.empty:
        mode_count = series.value_counts()[modes.iloc[0]]
        return mode_count
    else:
        return 0

# Calculate the average start and end positions for each group
result_df = df_IS.groupby('group').agg({'tseq': 'first', 'integration_variant':'first','start': get_mode, 
'end': get_mode,'Sample':'nunique','tstart_reads_evidence':'mean','tend_reads_evidence':'mean',"group":"first"}).reset_index(drop=True)
result_df['Prevalence']=result_df['Sample']/22
result_df.sort_values('Prevalence',ascending=False)  

# %%
samples_IS = {sample:[] for sample in df_IS['Sample'].unique()}
for group in result_df['group']:
    temp=df_IS[df_IS['group']==group]
    for sample in samples_IS.keys():
        if sample not in list(temp['Sample']):
            samples_IS[sample].append(0)
        else:
            samples_IS[sample].append(1)
for sample in samples_IS.keys():
    result_df[sample]=samples_IS[sample]
result_df

# %%
samples=list(all_integration_sites.columns)
for sample in samples:
    print(sample, ': ',str(len(result_df[(result_df[sample]==1)])))

# %%
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt

#distance_matrix = squareform(pdist(results_wild[samples_wild].T, metric='cityblock'))
distance_matrix = squareform(pdist(result_df[samples].T, metric='jaccard'))
linkage_matrix = linkage(distance_matrix, method='average')


# Custom dendrogram function with threshold line
def plot_dendrogram_with_threshold(linkage_matrix, threshold, ax=None):
    # Create dendrogram
    dendro = dendrogram(linkage_matrix, ax=ax, color_threshold=threshold,labels=samples,
    above_threshold_color='gray')
    
    # Add line at threshold
    plt.axhline(y=threshold, color='r', linestyle='--')
    plt.xticks(rotation=90, ha='right')
    # Color branches based on distance (auto-colored by dendrogram)
    # axhline has already added a red threshold line
    
    return dendro

# Plotting the custom dendrogram
fig, ax = plt.subplots(figsize=(15, 8))
dendro = plot_dendrogram_with_threshold(linkage_matrix, threshold=150, ax=ax)  # Set your threshold
ax.set_ylabel('Manhattan Distance')
#ax.set_xticklabels(samples_stbeespacbio)
plt.show()


# %%
#from sklearn.decomposition import PCA
pca = PCA(n_components=10)

x= pca.fit_transform(result_df[samples].T) ## PCA over samples
loadings = pd.DataFrame(x,
                        columns=['PC%s' % _ for _ in range(1,11)],
                        index=samples)
sns.lmplot(data=loadings,
    x="PC1", 
    y="PC2",  
    fit_reg=False, 
    legend=True)

    # Annotate each point with the sample name
for i, row in loadings.iterrows():
    plt.text(row['PC1'], row['PC2'], i, horizontalalignment='left', size='medium', color='black', weight='semibold')

plt.xlabel("PC1 ("+str(round(pca.explained_variance_ratio_[0]*100,2))+"%)")
plt.ylabel("PC2 ("+str(round(pca.explained_variance_ratio_[1]*100,2))+"%)")
plt.show()

# %%
import seaborn as sns
sns.boxplot(x="integration_variant", y="Prevalence", data=result_df)
sns.despine(offset=10, trim=True)

# %%
annot = pd.read_csv("GCA_030178435.1_masked_KoRV_PhER_PhaCinBeta/braker/braker.gtf",sep="\t",
names=['Contig','Source','Type','start','end','col','strand','col2','INFO'])
contigs={}
contigs_f = open("GCA_030178435.1/GCA_030178435.1_ASM3017843v1_genomic.fna.ann","r")
for line in contigs_f.readlines():
    if "whole genome shotgun sequence" in line:
        contig=line.split(" ")[6][:-1]
        name=line.split(" ")[1]
        contigs[contig]=name
annot['Contig']=annot['Contig'].replace(contigs)
contigs_f.close()
annot

# %%
orthologs = pd.read_excel("GCA_030178435.1/Supp_Table2.xlsx",comment="#")
orthologs = orthologs[~orthologs['query'].isna()]
orthologs['query']=orthologs['query'].apply(lambda x: x.split(".")[0])
genes = annot[annot['Type']=='gene'].reset_index(drop=True)
genes = pd.merge(orthologs,genes,how='left',left_on='query',right_on='INFO')
genes

# %%
overlapping_genes = []
flanking= 10000
genes['Preferred_name']=[genes.loc[i,'seed_ortholog'] if '-'==genes.loc[i,'Preferred_name'] else genes.loc[i,'Preferred_name'] for i in genes.index]
# Iterate over each row in contigs_df
for contig_row in result_df.itertuples(index=False):
    # Filter genes_df to find overlapping genes for the current contig
    overlapping_genes_contig = genes[
        (genes['Contig'] == contig_row.tseq) &
        (genes['start'] <= contig_row.end+flanking) &
        (genes['end'] >= contig_row.start-flanking)
    ]

    # Append the overlapping genes to the list
    overlapping_genes.append(', '.join(overlapping_genes_contig['Preferred_name'].tolist()))
result_df['Genes']=overlapping_genes
result_df

# %%
all_genes=[]
file = open("resultskorv_genes.txt","w")
for gene in result_df['Genes']:
    gene = gene.split(",")
    for item in gene:
        if item.strip() not in all_genes:
            all_genes.append(item.strip())
            file.write(item.strip()+"\n")
file.close()
len(all_genes)

# %%
all_genes=[]
for gene in result_df[result_df['Prevalence']>=0.4]['Genes']:
    gene = gene.split(",")
    for item in gene:
        if item.strip() not in all_genes:
            all_genes.append(item.strip())
            print(item.strip())

# %%
result_df[result_df['Genes'].str.contains('BCL2L1|SLC29A1|ALDH1A3|LZTS1')]

# %%
result_df[result_df['integration_variant']=='KoRV-B'].describe() #all KoRVB are somatic

# %%
df_IS['Number of reads'] = df_IS['tstart_reads_evidence'] + df_IS['tend_reads_evidence']
df_IS['AVG Number of reads'] = df_IS['Number of reads'] /2
df_IS['Max coverage'] = df_IS[['coverage_tstart','coverage_tend']].max(axis=1)
df_IS.describe()

# %%
df_IS['Region_coverage']=df_IS['Number of reads']/(df_IS['windown_len']+1)
df_IS['Region_coverage'].describe()

# %%
df_IS.columns

# %%
#get info per animal
sample_info = pd.read_excel('sample_info.xlsx')
sample_info['Number of IS']=[len(result_df[result_df[sample]==1]) for sample in sample_info['Name']]
sample_info['Number of KoRV-B']=[len(result_df[(result_df[sample]==1) & (result_df['integration_variant']=='KoRV-B')]) for sample in sample_info['Name']]
sample_info['Avg Number of reads']=[df_IS[df_IS['Sample']==sample]['Number of reads'].mean() for sample in sample_info['Name']]
sample_info['Stde number of reads']=[df_IS[df_IS['Sample']==sample]['Number of reads'].std() for sample in sample_info['Name']]
sample_info['Avg Max_cov']=[df_IS[df_IS['Sample']==sample]['Max coverage'].std() for sample in sample_info['Name']]
sample_info

# %%
sample_info.index=sample_info['Name']
df_IS=df_IS.reset_index(drop=True)
df_IS['Corrected_number_reads']=[df_IS.loc[i,'Number of reads']*(sample_info.loc[df_IS.loc[i,'Sample'],'Avg Number of reads']/sample_info['Avg Number of reads'].mean()) for i in df_IS.index]
df_IS['Corrected_max_coverage']=[df_IS.loc[i,'Max coverage']*(sample_info.loc[df_IS.loc[i,'Sample'],'Avg Max_cov']/sample_info['Avg Max_cov'].mean()) for i in df_IS.index]
df_IS

# %%
df_IS.index=df_IS["group"]
trios = {"koala triad names" }
sample_c=[];uniques=[];coverages=[]
for sample in trios.keys():
    temp=result_df[result_df[sample]>0]
    for i in temp.index:
        sample_c.append(sample)
        uniques.append(False if temp.loc[i,trios[sample][0]]==1 or  temp.loc[i,trios[sample][1]]==1 else True)
        coverages.append(int(df_IS[df_IS["Sample"]==sample].loc[temp.loc[i,"group"],"Corrected_number_reads"]))
df_covs=pd.DataFrame()
df_covs["Sample"]=sample_c
df_covs["Unique"]=uniques
df_covs["Coverage"]=coverages
df_covs.describe()

# %%
plt.figure(figsize=(10, 6))
g = sns.boxplot(x="Sample", y="Coverage",hue="Unique",
            data=df_covs)
plt.axhline(150, color='red', linestyle='--')
plt.yscale('log')
sns.despine(offset=10, trim=True)
plt.savefig('proportion_denovoVscommonIS.jpg', dpi=600)
plt.show()

# %%
sample_info.describe()

# %%
sample_info.to_excel('results_per_sample.xlsx',index=False)

# %%
df_IS.describe()

# %%
df_IS.groupby("Sample")["Number of reads"].sum()

# %%
sample_info=sample_info.astype({'Triad ':str})
for triad in [1,2,3,4,5,6]:#7 is not a triad, but just offspring to validate
    temp = sample_info[sample_info['Triad '].str.contains(str(triad))].reset_index(drop=True)
    temp['Triad function'] = [ x if '/' not in x else (x.split('/')[0] if str(triad)+'/' in temp.loc[i,'Triad '] else x.split('/')[1]) for i,x in enumerate(temp['Triad function'])] 
    dam = temp[temp['Triad function']=='dam']['Name'].iloc[0]
    sire = temp[temp['Triad function']=='sire']['Name'].iloc[0]
    for joey in temp[temp['Triad function']=='joey']['Name']:
        cov_threshold=500
        temp = result_df[(result_df[joey]==1) & (result_df[sire]==0) & (result_df[dam]==0) & (result_df['Sample']==1)]
        if joey=='395':
            temp = result_df[((result_df[joey]==1) & (result_df[sire]==0) & (result_df[dam]==0) & (result_df['Sample']==2) & (result_df['475']==1)) #IS present in F2 too
            | ((result_df[joey]==1) & (result_df[sire]==0) & (result_df[dam]==0) & (result_df['Sample']==1))] #IS present in F1 only
        elif joey=='369':
            cov_threshold=1000
            temp = result_df[((result_df[joey]==1) & (result_df[sire]==0) & (result_df[dam]==0) & (result_df['Sample']==2) & (result_df['530']==1)) |
            ((result_df[joey]==1) & (result_df[sire]==0) & (result_df[dam]==0) & (result_df['Sample']==1) )]
        elif joey=='241':
            cov_threshold=300
        temp['Coverage']=[df_IS[(df_IS['group']==group) & (df_IS['Sample']==joey)]['Number of reads'].iloc[0] for group in temp['group']]
        #cov_threshold = df_IS[(df_IS['Sample']==joey)]['Number of reads'].median()
        temp=temp[temp['Coverage']>=cov_threshold]
        print(str(triad)+' : '+str(len(temp)))
        print(temp[['tseq','start','group','Coverage','Genes',"Sample",joey,sire,dam]])

# %%
total_sequence_length = 2*(3234982288) #wilpena
5 / (total_sequence_length * 9)

# %%
total_sequence_length = 2*(3234982288) #wilpena
21 / (total_sequence_length * 9) #removing 373 and 369 - too many IS

# %%
5/9

# %%
joey='395';sire='201';dam='296' 
temp = result_df[(result_df[joey]==1) & (result_df[sire]==0) & (result_df[dam]==0) & (result_df['475']==1) & (result_df['Sample']==2)]
temp['Coverage']=[df_IS[(df_IS['group']==group) & (df_IS['Sample']==joey)]['Corrected_max_coverage'].iloc[0] for group in temp['group']]
temp=temp[temp['Coverage']>=(df_IS["Corrected_max_coverage"].mean())]
print(temp[['tseq','start','group','Coverage','Genes','Sample',joey,sire,dam,'475']])
