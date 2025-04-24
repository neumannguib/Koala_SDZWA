# %%
import glob
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="ticks", palette="colorblind")

# %%
#get sample info 
path = 'path to illumina data'
sample_info = pd.read_excel(path + 'Manifest.xlsx')
sample_info

# %%
result_df=pd.read_csv(path+"IS_detected_frequencies_annotated.tab",sep="\t")
result_df

# %%
samples_sandiego=[str(col) for col in sample_info["Sample ID"].unique()]
len(samples_sandiego)

# %%
pheno_SDZWA=pd.DataFrame()
pheno_SDZWA["Accession"]=samples_sandiego
temp=pd.read_csv(path+"IS_frequenciesv2.tab",sep="\t")
#consensus was considered for replicates
pheno_SDZWA['PhaCinBeta']=[len(temp[(temp['integration_variant']=='PhaCinBeta') & (temp[sample]>0)]) for sample in pheno_SDZWA['Accession']]
pheno_SDZWA['PhaCinBetalike']=[len(temp[(temp['integration_variant']=='PhaCinBetalike') & (temp[sample]>0)]) for sample in pheno_SDZWA['Accession']]

pheno_SDZWA['KoRV_filtered']=[len(result_df[(result_df['integration_variant']=='KoRV') & (result_df[sample]>0)]) for sample in pheno_SDZWA['Accession']]
pheno_SDZWA['PhaCinBeta_filtered']=[len(result_df[(result_df['integration_variant']=='PhaCinBeta') & (result_df[sample]>0)]) for sample in pheno_SDZWA['Accession']]
pheno_SDZWA['PhaCinBetalike_filtered']=[len(result_df[(result_df['integration_variant']=='PhaCinBetalike') & (result_df[sample]>0)]) for sample in pheno_SDZWA['Accession']]
pheno_SDZWA.describe()

# %%
results_wild = pd.read_csv(path+"IS_AAfrequencies_wild.tab",sep="\t")
results_wild

# %%
results_wild.groupby("integration_variant").count()['group']

# %%
results_EU=pd.read_excel(path+'IS_detected_annotatedEU.xlsx').replace({"KoRV-A":"KoRV","KoRV-B":"KoRV"})
results_EU

# %%
results_EU["Total reads"]=results_EU["tstart_reads_evidence"]+results_EU["tend_reads_evidence"]
results_EU[results_EU["Total reads"]>=100]

# %%
results_EU.describe()

# %%
samples_eu = list(results_EU.columns[-22:-2])
samples_eu 

# %%
len(samples_eu)

# %%
pheno_eu=pd.DataFrame()
pheno_eu["Accession"]=samples_eu
pheno_eu['KoRV']=[len(results_EU[(results_EU['integration_variant']=='KoRV') & (results_EU[str(sample)]>0)]) for sample in pheno_eu['Accession']]
pheno_eu['enKoRV']=[len(results_EU[(results_EU['integration_variant']=='KoRV') & (results_EU[str(sample)]>0) & (results_EU["Total reads"]>=100)]) for sample in pheno_eu['Accession']]
pheno_eu['Somatic KoRV']=[len(results_EU[(results_EU['integration_variant']=='KoRV') & (results_EU[str(sample)]>0) & (results_EU["Total reads"]<100)]) for sample in pheno_eu['Accession']]
pheno_eu.describe()

# %%
pheno_eu['KoRV_ERV']=[len(results_EU[(results_EU['integration_variant']=='KoRV') & (results_EU[str(sample)]>0) &
(results_EU["Total reads"]>=100) ]) for sample in pheno_eu['Accession']]
pheno_eu['KoRV_somatic']=[len(results_EU[(results_EU['integration_variant']=='KoRV') & (results_EU[str(sample)]>0) &
(results_EU["Total reads"]<100) ]) for sample in pheno_eu['Accession']]
pheno_eu.describe()

# %%
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
samples_wild=list(results_wild.columns[9:-6])
distance_matrix = squareform(pdist(results_wild[samples_wild].T, metric='jaccard'))
sns.clustermap(distance_matrix, cmap="viridis", method='average', figsize=(18, 18),
               yticklabels=samples_wild, xticklabels=samples_wild)

# %%
len(samples_wild)

# %%
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

#distance_matrix = squareform(pdist(results_wild[samples_wild].T, metric='cityblock'))
distance_matrix = squareform(pdist(results_wild[samples_wild].T, metric='jaccard'))
linkage_matrix = linkage(distance_matrix, method='average')


# Custom dendrogram function with threshold line
def plot_dendrogram_with_threshold(linkage_matrix, threshold, ax=None):
    # Create dendrogram
    dendro = dendrogram(linkage_matrix, ax=ax, color_threshold=threshold,labels=samples_wild,
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
distance_matrix = squareform(pdist(results_EU[samples_eu].T, metric='jaccard'))
linkage_matrix = linkage(distance_matrix, method='average')


# Custom dendrogram function with threshold line
def plot_dendrogram_with_threshold(linkage_matrix, threshold, ax=None):
    # Create dendrogram
    dendro = dendrogram(linkage_matrix, ax=ax, color_threshold=threshold,labels=samples_eu,
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
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

#distance_matrix = squareform(pdist(results_wild[samples_wild].T, metric='cityblock'))
distance_matrix = squareform(pdist(result_df[samples_sandiego].T, metric='jaccard'))
linkage_matrix = linkage(distance_matrix, method='average')


# Custom dendrogram function with threshold line
def plot_dendrogram_with_threshold(linkage_matrix, threshold, ax=None):
    # Create dendrogram
    dendro = dendrogram(linkage_matrix, ax=ax, color_threshold=threshold,labels=samples_sandiego,
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
starts = []; ends=[]
for i in results_EU.index:
    start=results_EU.loc[i,'start']
    end=results_EU.loc[i,'end']
    if start>end:
        temp=start
        start=end
        end=temp
    elif start==end:
        start=start-10
        end=end+10
    elif (end-start)<20:
        start=start-(20-(end-start))
        end=end+(20-(end-start))
    starts.append(start)
    ends.append(end)
results_EU['start']=starts
results_EU['end']=ends
results_EU['windown']=results_EU['end']-results_EU['start']
results_EU['windown'].describe()

# %%
import pandas as pd
from intervaltree import Interval, IntervalTree

# The three DataFrame sample lists
samples_all = samples_sandiego + samples_wild + samples_eu
global samples_all

# A function that processes overlaps
def find_overlaps(dfs):
    overlapped_rows = []
    overlap_indices = [set() for _ in range(len(dfs))]

    all_tseq_variant = set()
    for df in dfs:
        all_tseq_variant.update(tuple(row) for row in df[['tseq', 'integration_variant']].dropna().values)

    for (tseq, integration_variant) in all_tseq_variant:
        # Filter dataframes by tseq and integration_variant matches
        dfs_filtered = [df[(df['tseq'] == tseq) & (df['integration_variant'] == integration_variant)] for df in dfs]

        itrees = [IntervalTree() for _ in range(len(dfs_filtered))]
        #fill in a list of IS
        for i, df_filtered in enumerate(dfs_filtered):
            for idx, row in df_filtered.iterrows():
                start, end = row['start'], row['end']
                if pd.notnull(start) and pd.notnull(end):  # Check for null or NaN values in start/end columns
                    itrees[i].addi(start, end, idx)
        #check the overlaps
        for i, df_filtered in enumerate(dfs_filtered[:-1]):
            for idx, row in dfs_filtered[i].iterrows():
                overlaps_all=[]
                for j in range(i + 1, len(dfs_filtered)):
                    overlaps = itrees[j][row['start']:row['end']]
                    if overlaps:
                        for interval in overlaps:
                            overlap_data = {**dfs[i].loc[idx].to_dict(), **dfs[j].loc[interval.data].to_dict()}
                            # Handling different column names for start and end to avoid conflict
                            overlap_data.update({f'df{i+1}_start': row['start'], f'df{i+1}_end': row['end'], f'df{j+1}_start': interval.begin, f'df{j+1}_end': interval.end})
                            overlaps_all.append(overlap_data)
                            overlap_indices[j].add(interval.data)
                            overlap_indices[i].add(idx)
                if len(overlaps_all)==1:
                    overlapped_rows.append(overlaps_all[0])
                elif len(overlaps_all)>1:
                    temp=overlaps_all[0]
                    for k in range(1,len(overlaps_all)):
                        temp = {**temp,**overlaps_all[k]}
                    overlapped_rows.append(temp)
                

    # Make sure to rename start and end columns to prevent conflicts.
    unique_dfs = [df.loc[~df.index.isin(overlap_index)] for df, overlap_index in zip(dfs, overlap_indices)]
    for i, df in enumerate(unique_dfs):
        df.rename(columns={'start': f'df{i+1}_start', 'end': f'df{i+1}_end'}, inplace=True)

    # Creating a single merged DataFrame
    merged_df = pd.concat([pd.DataFrame(overlapped_rows)] + unique_dfs, ignore_index=True)

    return merged_df

# Merging the dataframes result_df, results_wild, and results_EU
results_EU['integration_variant']=results_EU['integration_variant'].replace({'KoRV-B':'KoRV','KoRV-A':'KoRV','PhaCinBetalike514':'PhaCinBetalike'})
merged_df = find_overlaps([result_df, results_wild, results_EU[(results_EU['integration_variant']=='KoRV') & (results_EU["Total reads"]>=100)].reset_index(drop=True)])
merged_df['Sample']=merged_df[samples_all].apply(lambda row: (row > 0).sum(), axis=1)
merged_df['start']=[merged_df.loc[i,'df1_start'] if merged_df['start'].isna().iloc[i] else merged_df.loc[i,'start'] for i in merged_df.index]
merged_df=merged_df.sort_values('Sample',ascending=False).drop_duplicates(['tseq','start','integration_variant','df2_start','df3_start'],keep='first') # remove duplicates caused by multiple overlaps
merged_df

# %%
merged_df.describe()

# %%
for col in merged_df.columns:
    print(col)

# %%
sample_info_wild = pd.read_excel("samples_selected.xlsx")
sample_info_wild.index=sample_info_wild["AWS File Name"]
sample_info_wild

# %%
merged_df[samples_wild]=merged_df[samples_wild].fillna(0)
merged_df[samples_eu]=merged_df[samples_eu].fillna(0)
merged_df[samples_sandiego]=merged_df[samples_sandiego].fillna(0)
samples_NSW =[sample for sample in samples_wild if sample_info_wild.loc[sample,"State"] =="NSW" ]
samples_QLD =[sample for sample in samples_wild if sample_info_wild.loc[sample,"State"] =="QLD" ]
samples_VIC =[sample for sample in samples_wild if sample_info_wild.loc[sample,"State"] =="VIC" ]
merged_df['AAF_SanDiego']=merged_df[samples_sandiego].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
merged_df['AAF_NSW']=merged_df[samples_NSW].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
merged_df['AAF_QLD']=merged_df[samples_QLD].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
merged_df['AAF_VIC']=merged_df[samples_VIC].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
merged_df['Missing_wild']= merged_df[samples_wild].apply(lambda row: (row < 0).sum()/75, axis=1)
merged_df.groupby("integration_variant").count()['group']

# %%
merged_df=merged_df[(merged_df["Missing_wild"]<=0.1) | (merged_df["AAF_SanDiego"]>0)]
merged_df.groupby("integration_variant").count()['group']

# %%
merged_df[(merged_df['AAF_SanDiego']>0) & (merged_df['AAF_QLD']>0) & (merged_df['AAF_NSW']>0) & (merged_df['AAF_VIC']>0) ]

# %%
len(samples_NSW)

# %%
pheno_NSW=pd.DataFrame()
pheno_NSW["Accession"]=samples_NSW
pheno_NSW['KoRV']=[len(results_wild[(results_wild['integration_variant']=='KoRV') & (results_wild[str(sample)]>0)]) for sample in pheno_NSW['Accession']]
pheno_NSW['PhaCinBeta']=[len(results_wild[(results_wild['integration_variant']=='PhaCinBeta') & (results_wild[str(sample)]>0)]) for sample in pheno_NSW['Accession']]
pheno_NSW['PhaCinBetalike']=[len(results_wild[(results_wild['integration_variant']=='PhaCinBetalike') & (results_wild[str(sample)]>0)]) for sample in pheno_NSW['Accession']]
count=[]; reference_length=3234982288 #wilpena
for sample in pheno_NSW["Accession"]:
    with open(path+'/mapped/'+str(sample)+'_sorted.bam_stats.txt', 'r') as file:
        for line in file:
            if "mapped" in line:  # Find total bases mapped
                total_mapped_bases = int(line.split(' ')[0])*146
                count.append(total_mapped_bases / reference_length)
                break
        
pheno_NSW['Coverage']=count
pheno_NSW.describe()

# %%
len(samples_QLD)

# %%
pheno_QLD=pd.DataFrame()
pheno_QLD["Accession"]=samples_QLD
pheno_QLD['KoRV']=[len(results_wild[(results_wild['integration_variant']=='KoRV') & (results_wild[str(sample)]>0)]) for sample in pheno_QLD['Accession']]
pheno_QLD['PhaCinBeta']=[len(results_wild[(results_wild['integration_variant']=='PhaCinBeta') & (results_wild[str(sample)]>0)]) for sample in pheno_QLD['Accession']]
pheno_QLD['PhaCinBetalike']=[len(results_wild[(results_wild['integration_variant']=='PhaCinBetalike') & (results_wild[str(sample)]>0)]) for sample in pheno_QLD['Accession']]
count=[]; reference_length=3234982288 #wilpena
for sample in pheno_QLD["Accession"]:
    with open(path+'mapped/'+str(sample)+'_sorted.bam_stats.txt', 'r') as file:
        for line in file:
            if "mapped" in line:  # Find total bases mapped
                total_mapped_bases = int(line.split(' ')[0])*146
                count.append(total_mapped_bases / reference_length)
                break
        
pheno_QLD['Coverage']=count
pheno_QLD.describe()

# %%
pheno_QLD.sort_values('PhaCinBeta')

# %%
len(samples_VIC)

# %%
pheno_VIC=pd.DataFrame()
pheno_VIC["Accession"]=samples_VIC
pheno_VIC['KoRV']=[len(results_wild[(results_wild['integration_variant']=='KoRV') & (results_wild[str(sample)]>0)]) for sample in pheno_VIC['Accession']]
pheno_VIC['PhaCinBeta']=[len(results_wild[(results_wild['integration_variant']=='PhaCinBeta') & (results_wild[str(sample)]>0)]) for sample in pheno_VIC['Accession']]
pheno_VIC['PhaCinBetalike']=[len(results_wild[(results_wild['integration_variant']=='PhaCinBetalike') & (results_wild[str(sample)]>0)]) for sample in pheno_VIC['Accession']]
count=[]; reference_length=3234982288 #wilpena
for sample in pheno_VIC["Accession"]:
    with open(path+'mapped/'+str(sample)+'_sorted.bam_stats.txt', 'r') as file:
        for line in file:
            if "mapped" in line:  # Find total bases mapped
                total_mapped_bases = int(line.split(' ')[0])*146
                count.append(total_mapped_bases / reference_length)
                break
        
pheno_VIC['Coverage']=count
pheno_VIC.describe()

# %%
pheno_VIC[["KoRV","PhaCinBeta", "PhaCinBetalike", "Coverage"]].corr()

# %%
from scipy.stats import pearsonr
correlation_matrix = pheno_VIC[["KoRV","PhaCinBeta", "PhaCinBetalike", "Coverage"]].corr()
columns =["KoRV","PhaCinBeta", "PhaCinBetalike", "Coverage"]
# Initialize a DataFrame to store p-values
p_values = pd.DataFrame(np.ones(correlation_matrix.shape), columns=columns, index=columns)

# Compute p-values
for col1 in columns:
    for col2 in columns:
        if col1 != col2:
            corr, p_val = pearsonr(pheno_VIC[col1], pheno_VIC[col2])
            p_values.at[col1, col2] = p_val
        else:
            p_values.at[col1, col2] = np.nan  # NaN for self-correlation
p_values

# %%
244/1232

# %%
pheno_NSW[["KoRV","PhaCinBeta", "PhaCinBetalike", "Coverage"]].corr()

# %%
pheno_NSW

# %%
correlation_matrix = pheno_NSW[["KoRV","PhaCinBeta", "PhaCinBetalike", "Coverage"]].corr()
columns =["KoRV","PhaCinBeta", "PhaCinBetalike", "Coverage"]
# Initialize a DataFrame to store p-values
p_values = pd.DataFrame(np.ones(correlation_matrix.shape), columns=columns, index=columns)

# Compute p-values
for col1 in columns:
    for col2 in columns:
        if col1 != col2:
            corr, p_val = pearsonr(pheno_NSW[col1], pheno_NSW[col2])
            p_values.at[col1, col2] = p_val
        else:
            p_values.at[col1, col2] = np.nan  # NaN for self-correlation
p_values

# %%
pheno_QLD[["KoRV","PhaCinBeta", "PhaCinBetalike", "Coverage"]].corr()

# %%
from scipy.stats import pearsonr
correlation_matrix = pheno_QLD[["KoRV","PhaCinBeta", "PhaCinBetalike", "Coverage"]].corr()
columns =["KoRV","PhaCinBeta", "PhaCinBetalike", "Coverage"]
# Initialize a DataFrame to store p-values
p_values = pd.DataFrame(np.ones(correlation_matrix.shape), columns=columns, index=columns)

# Compute p-values
for col1 in columns:
    for col2 in columns:
        if col1 != col2:
            corr, p_val = pearsonr(pheno_QLD[col1], pheno_QLD[col2])
            p_values.at[col1, col2] = p_val
        else:
            p_values.at[col1, col2] = np.nan  # NaN for self-correlation
p_values

# %%
samples_minusEU = [x for x in samples_all if x not in samples_eu]
merged_df['Missing']= merged_df[samples_minusEU].apply(lambda row: (row < 0).sum()/len(samples_minusEU), axis=1)
merged_df[merged_df["Missing"]<=0.1].groupby("integration_variant").count()['group']

# %%
merged_df[(merged_df['AAF_SanDiego']>0) & (merged_df['AAF_QLD']>0) & (merged_df['AAF_NSW']>0) & (merged_df['AAF_VIC']>0)].groupby("integration_variant").count()["group"]

# %%
merged_df[(merged_df['AAF_SanDiego']>0) & (merged_df['AAF_QLD']>0) & (merged_df['AAF_NSW']>0) 
& (merged_df['AAF_VIC']>0) & (merged_df["integration_variant"]=="PhaCinBeta")][["tseq","start","end","POS","AAF_SanDiego","AAF_NSW","AAF_QLD","AAF_VIC","Genes"]]

# %%
mask = (merged_df[samples_all]>=0).all(axis=1)
merged_df_nomissing = merged_df[mask]
merged_df['NSW (n=25)']=merged_df[samples_NSW].apply(lambda row: (row > 0).sum(), axis=1)
merged_df['QLD (n=25)']=merged_df[samples_QLD].apply(lambda row:  (row > 0).sum(), axis=1)
merged_df['VIC (n=25)']=merged_df[samples_VIC].apply(lambda row:  (row > 0).sum(), axis=1)
merged_df['SDZWA (n=91)']=merged_df[samples_sandiego].apply(lambda row:  (row > 0).sum(), axis=1)
merged_df['EUZ (n=20)']=merged_df[samples_eu].apply(lambda row:  (row > 0).sum(), axis=1)
merged_df=merged_df[(merged_df['NSW (n=25)']>0) | (merged_df['QLD (n=25)']>0) | (merged_df['VIC (n=25)']>0) | (merged_df['SDZWA (n=91)']>0)
| (merged_df['EUZ (n=20)']>0)].reset_index(drop=True)

merged_df_nomissing = merged_df[(merged_df['AAF_SanDiego']>=0.05) | (merged_df['AAF_QLD']>=0.05) | (merged_df['AAF_NSW']>=0.05) 
| (merged_df['AAF_VIC']>=0.05)][samples_all].fillna(0).replace({-2:0,-1:0,-3:0,-4:0,-5:0})

merged_df_nomissing

# %%
merged_df[(merged_df['AAF_SanDiego']>0)].groupby("integration_variant").count()["group"]

# %%
merged_df[(merged_df['EUZ (n=20)']>0)].groupby("integration_variant").count()["group"]

# %%
merged_df[(merged_df['AAF_QLD']>0) | (merged_df['AAF_NSW']>0) | (merged_df['AAF_VIC']>0)].groupby("integration_variant").count()["group"]

# %%
len(merged_df[(merged_df['AAF_QLD']>0) | (merged_df['AAF_NSW']>0) | (merged_df['AAF_VIC']>0)])

# %%
temp=merged_df[(merged_df['AAF_NSW']>0) ]
pheno_NSW['KoRV_filtered']=[len(temp[(temp['integration_variant']=='KoRV') & (temp[str(sample)]>0)]) for sample in pheno_NSW['Accession']]
pheno_NSW['PhaCinBeta_filtered']=[len(temp[(temp['integration_variant']=='PhaCinBeta') & (temp[str(sample)]>0)]) for sample in pheno_NSW['Accession']]
pheno_NSW['PhaCinBetalike_filtered']=[len(temp[(temp['integration_variant']=='PhaCinBetalike') & (temp[str(sample)]>0)]) for sample in pheno_NSW['Accession']]
pheno_NSW.describe()

# %%
temp=merged_df[(merged_df['AAF_QLD']>0) ]
pheno_QLD['KoRV_filtered']=[len(temp[(temp['integration_variant']=='KoRV') & (temp[str(sample)]>0)]) for sample in pheno_QLD['Accession']]
pheno_QLD['PhaCinBeta_filtered']=[len(temp[(temp['integration_variant']=='PhaCinBeta') & (temp[str(sample)]>0)]) for sample in pheno_QLD['Accession']]
pheno_QLD['PhaCinBetalike_filtered']=[len(temp[(temp['integration_variant']=='PhaCinBetalike') & (temp[str(sample)]>0)]) for sample in pheno_QLD['Accession']]
pheno_QLD.describe()

# %%
temp=merged_df[(merged_df['AAF_VIC']>0) ]
pheno_VIC['KoRV_filtered']=[len(temp[(temp['integration_variant']=='KoRV') & (temp[str(sample)]>0)]) for sample in pheno_VIC['Accession']]
pheno_VIC['PhaCinBeta_filtered']=[len(temp[(temp['integration_variant']=='PhaCinBeta') & (temp[str(sample)]>0)]) for sample in pheno_VIC['Accession']]
pheno_VIC['PhaCinBetalike_filtered']=[len(temp[(temp['integration_variant']=='PhaCinBetalike') & (temp[str(sample)]>0)]) for sample in pheno_VIC['Accession']]
pheno_VIC.describe()

# %%
merged_df.groupby("integration_variant").count()['group']

# %%
merged_df.to_csv(path+"merged_IS_SanDiego_wild.tab",index=False,sep="\t")

# %%
#save as vcf
merged_df['FILTER'] = 'PASS'
merged_df['QUAL'] = 40
merged_df['FORMAT']='GT'
merged_df['INFO']="."
merged_df['is_snp']=False
merged_df['ALT'] = '<INS>'
merged_df['REF'] ='N'
merged_df['CHROMPOS']=merged_df["tseq"] + '_' + merged_df["POS"].astype(str) + merged_df["integration_variant"]
geno= merged_df[['tseq', 'POS', 'CHROMPOS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT','is_snp']+samples_all] 
geno=geno[~geno['POS'].isna()]
geno['Position']=geno['POS'].apply(lambda x: int(x) if '-' not in x else int(x.split('-')[0]))
geno=geno.sort_values(by=['tseq', 'Position'])
header = """##fileformat=VCFv4.2
##fileDate=05.09.2024
##source=retroviral insertions
##reference=GCA_030178435.1
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=Variant_type,Number=1,Type=String,Description="dbSNP ID">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"""

for col in samples_all:
    geno[col]=geno[col].replace({0:'0/0',1:'0/1',2:'1/1',-1:'./.',None:'./.',np.nan:'./.',-2:'./.'})
    header=header+'\t'+col
output_VCF = path+"merged_IS_SanDiego_wild.vcf"
geno=geno.reset_index(drop=True)

with open(output_VCF, 'w') as vcf:
    vcf.write(header)
    vcf.write('\n')
geno[['tseq', 'POS', 'CHROMPOS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']+samples_all].to_csv(output_VCF, sep="\t", mode='a', index=False,header=False)
vcf.close()

# %%
os.system("bgzip path_to_merged_IS_SanDiego_wild.vcf -f")
os.system("tabix path_to_merged_IS_SanDiego_wild.vcf.gz")

# %%
len(merged_df)

# %%
merged_df[merged_df['integration_variant']=='KoRV'].sort_values('AAF',ascending=False)[['AAF','AAF-alive','AAF-dead','AAF-diff',
'AAF_SanDiego','AAF_NSW','AAF_QLD','AAF_VIC','tseq','start','end','POS','Genes','EUZ (n=20)','df2_start','df2_end','df3_start','df3_end']]

# %%
merged_df[(merged_df['integration_variant']=='KoRV') & (~merged_df['Genes'].isna()) & (merged_df['Sample']>20)].sort_values('AAF',ascending=False)[['AAF','AAF-alive','AAF-dead','AAF-diff',
'AAF_SanDiego','AAF_NSW','AAF_QLD','AAF_VIC','tseq','start','end','POS','Genes','EUZ (n=20)','df2_start','df2_end','df3_start','df3_end']]

# %%
#vep --input_file path_to_merged_IS_SanDiego_wild.vcf.gz --output_file merged_IS_SanDiego_wild.annotated.vcf --fasta GCA_030178435.1_ASM3017843v1_genomic.fna --gff sorted_genomic_biotype.gff.gz -vcf 

# %%
def read_vcf(file_path):
    # Read the VCF file, skipping lines starting with '#'
    with open(file_path) as vcf_file:
        # Use the header to get column names from the VCF
        header = [line.strip() for line in vcf_file if line.startswith("#CHROM")]
        column_names =  header[0].split("\t")

    # Load the data into a DataFrame
    vcf_data = pd.read_csv(file_path, comment='#', sep='\t', names=column_names)

    return vcf_data

vcf_file_path = path+'merged_IS_SanDiego_wild.annotated.vcf'  
vcf_df = read_vcf(vcf_file_path)
vcf_df['Gene']=vcf_df["INFO"].apply(lambda x: x.split("|")[3])
vcf_df['Impact']=vcf_df["INFO"].apply(lambda x: x.split("|")[2])
vcf_df['Consequence']=vcf_df["INFO"].apply(lambda x: x.split("|")[1])
vcf_df.head()

# %%
vcf_df=vcf_df.drop_duplicates('ID',keep='first')
merged_df=pd.merge(merged_df,vcf_df[['ID','INFO','Impact','Consequence','Gene']],how='left',left_on='CHROMPOS',right_on='ID')
merged_df.head()

# %%
genes=merged_df['Gene'].unique()
len(genes)

# %%
merged_df['Consequence'].unique()

# %%
merged_df[merged_df['Consequence']=='feature_elongation']

# %%
merged_df['Consequence']=merged_df['Consequence'].fillna('')
merged_df[merged_df['Consequence'].str.contains("coding_sequence_variant")]

# %%
merged_df[merged_df['Consequence'].str.contains("intron_variant")]

# %%
len(merged_df)-1273+27

# %%
merged_df[(merged_df['Impact']=="HIGH") | (merged_df['Impact']=="MODERATE")][["tseq","POS","AAF_Total","AAF_SanDiego","EUZ (n=20)",
"SDZWA (n=91)","Genes","VIC (n=25)","QLD (n=25)","NSW (n=25)","Consequence"]]

# %%
for i,row in enumerate(merged_df.loc[3313,samples_sandiego]):
    if row==1:
        print(samples_sandiego[i])

# %%
for i,row in enumerate(merged_df.loc[810,samples_all]):#MCM9
    if row==1:
        print(samples_all[i])
    elif row==2:
        print(samples_all[i]+" HOMO")

# %%
merged_df[(merged_df['Impact']=="LOW") ]

# %%
merged_df[(merged_df['Gene']=="BCL2L1") ]

# %%
merged_df.loc[1517,"INFO_y"]

# %%
merged_df['Consequence']=merged_df['Consequence'].fillna('')
merged_df[(merged_df['Consequence'].str.contains("feature_elongation")) ]

# %%
merged_df[merged_df['Gene']=='LONP2'][["tseq","POS","AAF_Total","AAF_SanDiego","EUZ (n=20)",
"SDZWA (n=91)","Genes","VIC (n=25)","QLD (n=25)","NSW (n=25)"]]

# %%
merged_df[merged_df['Gene']=='STAP1'][["tseq","POS","AAF_Total","AAF_SanDiego","EUZ (n=20)","SDZWA (n=91)","Genes","VIC (n=25)","QLD (n=25)","NSW (n=25)"]]

# %%
import matplotlib.colors as mcolors
distance_matrix = squareform(pdist(merged_df_nomissing[samples_minusEU].T, metric='cityblock'))
linkage_matrix = linkage(distance_matrix, method='average')


# Plotting the custom dendrogram
fig, ax = plt.subplots(figsize=(7, 25))
dendro = dendrogram(linkage_matrix, ax=ax, labels=samples_minusEU,color_threshold=0, above_threshold_color="black",orientation='right')

def get_label_color(label_text):
    if label_text in list(samples_eu):
        return
    elif label_text in list(samples_sandiego):
        return mcolors.to_rgba("mediumseagreen")
    elif sample_info_wild.loc[label_text,"State"]=="NSW":
        return mcolors.to_rgba("orange")
    elif sample_info_wild.loc[label_text,"State"]=="QLD":
        return mcolors.to_rgba("lightblue")
    else:
        return mcolors.to_rgba("mistyrose")

        
# Set background color for each label in the dendrogram
for label in plt.gca().get_yticklabels():
    label.set_backgroundcolor(get_label_color(label.get_text()))
    label.set_fontsize(5)
    #label.set_bbox({'facecolor': get_label_color(label.get_text()), 'edgecolor': 'black'})

ax.set_xlabel('Manhattan Distance')
plt.tight_layout()  # Adjust layout for consistent background sizes

plt.show()

# %%
from sklearn.decomposition import PCA
plt.rcParams.update({
    'axes.titlesize': 14,     # Plot title font size
    'axes.labelsize': 20,     # X and Y axis labels font size
    'xtick.labelsize': 20,    # X tick labels font size
    'ytick.labelsize': 20,    # Y tick labels font size
    'legend.fontsize': 20,    # Legend font size
    'font.size': 20           # Base font size for other text
})
pca = PCA(n_components=10)
x= pca.fit_transform(merged_df_nomissing[samples_minusEU].T) ## PCA over samples
loadings = pd.DataFrame(x,
                        columns=['PC%s' % _ for _ in range(1,11)],
                        index=list(samples_minusEU))
plt.figure(figsize=(12,6))
plt.plot(pca.explained_variance_ratio_)
plt.ylabel('Explained Variance')
plt.xlabel('Components')
plt.xticks(list(range(10)),list(range(1,11)))
plt.show()

# %%
## PLOT PCA
pop = []
for sample in loadings.index:
    if sample in list(samples_sandiego):
        pop.append("SDZWA")
    else:
        pop.append(sample_info_wild.loc[sample,"State"])
loadings["Origin"]=pop
plt.figure(figsize=(10, 10))
sns.lmplot(
    x="PC1", 
    y="PC2", 
    data=loadings, 
    hue='Origin', 
    fit_reg=False, 
    legend=True,scatter_kws={'alpha': 1,'s': 100, 'edgecolor': 'black'},
    palette={"SDZWA":"mediumseagreen","QLD":"lightblue","NSW":"orange","VIC":"mistyrose"})

plt.xlabel("PC1 ("+str(round(pca.explained_variance_ratio_[0]*100,2))+"%)")
plt.ylabel("PC2 ("+str(round(pca.explained_variance_ratio_[1]*100,2))+"%)")
#plt.savefig('/home/guilherme/analyses/San_Diego/results/PCA_SanDiego_wildv2.jpg', dpi=1200)
plt.show()


# %%
sns.lmplot(
    x="PC1", 
    y="PC2", 
    data=loadings, 
    hue='Origin', 
    fit_reg=False, 
    legend=True,scatter_kws={'alpha': 1,'s': 100, 'edgecolor': 'black'},
    palette={"SDZWA":"mediumseagreen","QLD":"lightblue","NSW":"orange","VIC":"mistyrose"})
for i in range(loadings.shape[0]):
    if loadings["Origin"].iloc[i]=="SDZWA":
        plt.text(
            loadings["PC1"].iloc[i],
            loadings["PC2"].iloc[i],
            loadings.index[i],  # or use another column for labels
            fontsize=9,  # Adjust font size if necessary
            ha='right'   # Horizontal alignment, adjust as needed
        )
plt.xlabel("PC1 ("+str(round(pca.explained_variance_ratio_[0]*100,2))+"%)")
plt.ylabel("PC2 ("+str(round(pca.explained_variance_ratio_[1]*100,2))+"%)")
plt.show()

# %%
loadings

# %%
from sklearn.preprocessing import StandardScaler
pca = PCA(n_components=10)
merged_df["Total_Count"]=merged_df["NSW (n=25)"]+merged_df["QLD (n=25)"]+merged_df["VIC (n=25)"]+merged_df['SDZWA (n=91)']+merged_df['EUZ (n=20)']

merged_df_nomissing = merged_df[(merged_df['Total_Count']>=2) & (merged_df["integration_variant"]=="KoRV")][samples_all].fillna(0).replace({-2:0,-1:0,-3:0,-4:0,-5:0,2:1})
scaler = StandardScaler()
scaled_data=scaler.fit_transform(merged_df_nomissing[samples_all].T)

# Fit and transform the data
x= pca.fit_transform(scaled_data) ## PCA over samples

loadings = pd.DataFrame(x,
                        columns=['PC%s' % _ for _ in range(1,11)],
                        index=list(samples_all))
pop2=[]
for sample in loadings.index:
    if sample in list(samples_sandiego):
        pop2.append("SDZWA")
    elif sample in list(samples_eu):
        pop2.append("EUZ")
    else:
        pop2.append(sample_info_wild.loc[sample,"State"])
    
loadings["Origin"]=pop2
plt.figure(figsize=(10, 10))
sns.lmplot(
    x="PC1", 
    y="PC2", 
    data=loadings, 
    hue='Origin', 
    fit_reg=False, 
    legend=True,scatter_kws={'alpha': 1,'s': 100, 'edgecolor': 'black'},
    palette={"SDZWA":"mediumseagreen","QLD":"lightblue","NSW":"orange","VIC":"mistyrose","EUZ":"gray"})
plt.xlabel("PC1 ("+str(round(pca.explained_variance_ratio_[0]*100,2))+"%)")
plt.ylabel("PC2 ("+str(round(pca.explained_variance_ratio_[1]*100,2))+"%)")
plt.savefig('PCA_SanDiego_wild_EUZoos.jpg', dpi=1200)
plt.show()

# %%
len(merged_df[(merged_df['Total_Count']>=2) & (merged_df["integration_variant"]=="KoRV")])

# %%
plt.figure(figsize=(10, 10))
sns.lmplot(
    x="PC3", 
    y="PC4", 
    data=loadings, 
    hue='Origin', 
    fit_reg=False, 
    legend=True,scatter_kws={'alpha': 1,'s': 100, 'edgecolor': 'black'},
    palette={"SDZWA":"mediumseagreen","QLD":"lightblue","NSW":"orange","VIC":"mistyrose","EUZ":"gray"})
plt.xlabel("PC3 ("+str(round(pca.explained_variance_ratio_[2]*100,2))+"%)")
plt.ylabel("PC4 ("+str(round(pca.explained_variance_ratio_[3]*100,2))+"%)")
plt.savefig('PCA_SanDiego_wild_EUZoos_PC3AND4v2.jpg', dpi=1200)
plt.show()

# %%
plt.figure(figsize=(12,6))
plt.plot(pca.explained_variance_ratio_)
plt.ylabel('Explained Variance')
plt.xlabel('Components')
plt.xticks(list(range(10)),list(range(1,11)))
plt.show()

# %%
len(merged_df_nomissing)

# %%
pca = PCA(n_components=10)
merged_df_nomissing = merged_df[((merged_df['AAF_SanDiego']>=0.05) | (merged_df['AAF_QLD']>=0.05) | (merged_df['AAF_NSW']>=0.05) 
| (merged_df['AAF_VIC']>=0.05)) & (merged_df["integration_variant"]=="PhaCinBeta")][samples_minusEU].fillna(0).replace({-2:0,-1:0,-3:0,-4:0,-5:0})
x= pca.fit_transform(merged_df_nomissing[samples_minusEU].T) ## PCA over samples
loadings = pd.DataFrame(x,
                        columns=['PC%s' % _ for _ in range(1,11)],
                        index=list(samples_minusEU))
loadings["Origin"]=pop
sns.lmplot(
    x="PC1", 
    y="PC2", 
    data=loadings, 
    hue='Origin', 
    fit_reg=False, 
    legend=True)
plt.xlabel("PC1 ("+str(round(pca.explained_variance_ratio_[0]*100,2))+"%)")
plt.ylabel("PC2 ("+str(round(pca.explained_variance_ratio_[1]*100,2))+"%)")
plt.show()

# %%
pca = PCA(n_components=10)
merged_df_nomissing = merged_df[((merged_df['AAF_SanDiego']>=0.05) | (merged_df['AAF_QLD']>=0.05) | (merged_df['AAF_NSW']>=0.05) 
| (merged_df['AAF_VIC']>=0.05)) & (merged_df["integration_variant"]=="PhaCinBetalike")][samples_minusEU].fillna(0).replace({-2:0,-1:0,-3:0,-4:0,-5:0})
x= pca.fit_transform(merged_df_nomissing[samples_minusEU].T) ## PCA over samples
loadings = pd.DataFrame(x,
                        columns=['PC%s' % _ for _ in range(1,11)],
                        index=list(samples_minusEU))
loadings["Origin"]=pop
sns.lmplot(
    x="PC1", 
    y="PC2", 
    data=loadings, 
    hue='Origin', 
    fit_reg=False, 
    legend=True)
plt.xlabel("PC1 ("+str(round(pca.explained_variance_ratio_[0]*100,1))+"%)")
plt.ylabel("PC2 ("+str(round(pca.explained_variance_ratio_[1]*100,1))+"%)")
plt.show()

# %%
pca = PCA(n_components=10)
merged_df_nomissing = merged_df[((merged_df['AAF_SanDiego']>=0.05) | (merged_df['AAF_QLD']>=0.05) | (merged_df['AAF_NSW']>=0.05) 
| (merged_df['AAF_VIC']>=0.05))][samples_minusEU].fillna(0).replace({-2:0,-1:0,-3:0,-4:0,-5:0})
x= pca.fit_transform(merged_df_nomissing[samples_minusEU].T) ## PCA over samples
loadings = pd.DataFrame(x,
                        columns=['PC%s' % _ for _ in range(1,11)],
                        index=list(samples_minusEU))
loadings["Origin"]=pop
sns.lmplot(
    x="PC1", 
    y="PC2", 
    data=loadings, 
    hue='Origin', 
    fit_reg=False, 
    legend=True)
plt.xlabel("PC1 ("+str(round(pca.explained_variance_ratio_[0]*100,2))+"%)")
plt.ylabel("PC2 ("+str(round(pca.explained_variance_ratio_[1]*100,2))+"%)")
plt.show()

# %%
pca = PCA(n_components=10)
merged_df_nomissing = merged_df[( (merged_df['AAF_QLD']>=0.05) | (merged_df['AAF_NSW']>=0.05) 
| (merged_df['AAF_VIC']>=0.05))][samples_wild].fillna(0).replace({-2:0,-1:0,-3:0,-4:0,-5:0})
x= pca.fit_transform(merged_df_nomissing[samples_wild].T) ## PCA over samples
loadings = pd.DataFrame(x,
                        columns=['PC%s' % _ for _ in range(1,11)],
                        index=list(samples_wild))
pop = []
for sample in loadings.index:
    pop.append(sample_info_wild.loc[sample,"State"])

loadings["Origin"]=pop
sns.lmplot(
    x="PC1", 
    y="PC2", 
    data=loadings, 
    hue='Origin', 
    fit_reg=False, 
    legend=True)
plt.xlabel("PC1 ("+str(round(pca.explained_variance_ratio_[0]*100,2))+"%)")
plt.ylabel("PC2 ("+str(round(pca.explained_variance_ratio_[1]*100,2))+"%)")
plt.show()

# %%
pca = PCA(n_components=10)
merged_df_nomissing = merged_df[( (merged_df['AAF_SanDiego']>=0.05) | (merged_df['EUZ (n=20)']>=1) 
) & (merged_df["integration_variant"]=="KoRV")][samples_eu + samples_sandiego].fillna(0).replace({-2:0,-1:0,-3:0,-4:0,-5:0,2:1})
x= pca.fit_transform(merged_df_nomissing[samples_eu + samples_sandiego].T) ## PCA over samples
loadings = pd.DataFrame(x,
                        columns=['PC%s' % _ for _ in range(1,11)],
                        index=samples_eu + samples_sandiego)
pop = []
for sample in loadings.index:
    if sample in list(samples_eu):
        pop.append("EUZ")
    else:
        pop.append("SDZwa")

loadings["Origin"]=pop
sns.lmplot(
    x="PC1", 
    y="PC2", 
    data=loadings, 
    hue='Origin', 
    fit_reg=False, 
    legend=True)
plt.xlabel("PC1 ("+str(round(pca.explained_variance_ratio_[0]*100,2))+"%)")
plt.ylabel("PC2 ("+str(round(pca.explained_variance_ratio_[1]*100,2))+"%)")
plt.show()

# %%
#PLOTS - UNIQUE VERSUS SHARED, AND DISTRIBUTIONS
import pandas as pd
from upsetplot import UpSet
from upsetplot import plot
from upsetplot import generate_counts
plt.rcdefaults()
temp=merged_df[(merged_df['NSW (n=25)']>0) | (merged_df['QLD (n=25)']>0) | (merged_df['VIC (n=25)']>0) | (merged_df['SDZWA (n=91)']>0)]
temp['AAF']=temp[samples_minusEU].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
df_bool = temp[['QLD (n=25)','NSW (n=25)','VIC (n=25)','SDZWA (n=91)']] > 0
df_bool=df_bool.set_index(list(df_bool.columns))
df_bool['Insertion type']=list(temp['integration_variant'].replace({'PhaCinBeta':r'phaCin-$\beta$','PhaCinBetalike':r'phaCin-$\beta$-like'}))
df_bool['Total AAF']=list(temp['AAF'])
df_bool 
plt.figure(figsize=(50, 10))
upset = UpSet(df_bool, intersection_plot_elements=0,show_counts=True, shading_color=0.01,element_size=45)
upset.add_stacked_bars(by="Insertion type", colors=['indianred',sns.color_palette("husl", 8)[1],'saddlebrown'], title="Intersection Size", elements=4)
upset.add_catplot(value="Total AAF", kind="box", color="darkgray")
upset.style_categories('QLD (n=25)', bar_facecolor="lightblue",  bar_edgecolor="black")
upset.style_categories('NSW (n=25)', bar_facecolor="orange",  bar_edgecolor="black")
upset.style_categories('VIC (n=25)', bar_facecolor="mistyrose",bar_edgecolor="black")
upset.style_categories('SDZWA (n=91)', bar_facecolor="mediumseagreen", bar_edgecolor="black")
upset.plot()
plt.title("Intersection for KoRV, "+r'phaCin-$\beta$'+", and "+r'phaCin-$\beta$-like'+" Insertions")
plt.savefig('upsetplot_allISv2.jpg', dpi=1200)

plt.show()

# %%
(140+61+22+2)/2075

# %%
df_bool = merged_df[merged_df['integration_variant']=='KoRV'][['QLD (n=25)','NSW (n=25)','VIC (n=25)','SDZWA (n=91)','EUZ (n=20)']] > 0
df_bool=df_bool.set_index(list(df_bool.columns))
df_bool["Insertion type"]='KoRV'
df_bool 
plt.figure(figsize=(50, 10))
upset = UpSet(df_bool, intersection_plot_elements=0,show_counts=True, shading_color=0.01,element_size=45)
upset.add_stacked_bars(by="Insertion type", colors=['indianred'], title="Intersection Size", elements=5)


upset.style_categories('QLD (n=25)', bar_facecolor="lightblue",  bar_edgecolor="black")
upset.style_categories('NSW (n=25)', bar_facecolor="orange",  bar_edgecolor="black")
upset.style_categories('VIC (n=25)', bar_facecolor="mistyrose",bar_edgecolor="black")
upset.style_categories('SDZWA (n=91)', bar_facecolor="mediumseagreen", bar_edgecolor="black")
upset.style_categories('EUZ (n=20)', bar_facecolor="gray", bar_edgecolor="black")


upset.plot()
plt.title("Intersection for KoRV Insertions Only")
plt.savefig('upsetplot_KoRVv2.jpg', dpi=1200)

plt.show()

# %%
df_bool = merged_df[merged_df['integration_variant']=='PhaCinBeta'][['QLD (n=25)','NSW (n=25)','VIC (n=25)','SDZWA (n=91)']] > 0
df_bool=df_bool.set_index(list(df_bool.columns))
df_bool["Insertion type"]=r'phaCin-$\beta$'

plt.figure(figsize=(50, 10))
upset = UpSet(df_bool , intersection_plot_elements=0,show_counts=True, shading_color=0.01,element_size=45)
upset.add_stacked_bars(by="Insertion type", colors=[sns.color_palette("husl", 8)[1]], title="Intersection Size", elements=4)

upset.style_categories('QLD (n=25)', bar_facecolor="lightblue",  bar_edgecolor="black")
upset.style_categories('NSW (n=25)', bar_facecolor="orange",  bar_edgecolor="black")
upset.style_categories('VIC (n=25)', bar_facecolor="mistyrose",bar_edgecolor="black")
upset.style_categories('SDZWA (n=91)', bar_facecolor="mediumseagreen", bar_edgecolor="black")

upset.plot()
plt.title("Intersection for "+r'phaCin-$\beta$'+' Only')
plt.savefig('upsetplot_phacinbetav2.jpg', dpi=1200)

plt.show()

# %%
df_bool = merged_df[merged_df['integration_variant']=='PhaCinBetalike'][['QLD (n=25)','NSW (n=25)','VIC (n=25)','SDZWA (n=91)']] > 0
df_bool=df_bool.set_index(list(df_bool.columns))
df_bool["Insertion type"]=r'phaCin-$\beta$'+'-like'

plt.figure(figsize=(50, 10))
upset = UpSet(df_bool , intersection_plot_elements=0,show_counts=True, shading_color=0.01,element_size=45)
upset.add_stacked_bars(by="Insertion type", colors=['saddlebrown'], title="Intersection Size", elements=4)

upset.style_categories('QLD (n=25)', bar_facecolor="lightblue",  bar_edgecolor="black")
upset.style_categories('NSW (n=25)', bar_facecolor="orange",  bar_edgecolor="black")
upset.style_categories('VIC (n=25)', bar_facecolor="mistyrose",bar_edgecolor="black")
upset.style_categories('SDZWA (n=91)', bar_facecolor="mediumseagreen", bar_edgecolor="black")

upset.plot()
plt.title("Intersection for "+r'phaCin-$\beta$'+'-like Only')
plt.savefig('upsetplot_phacinbetalikev2.jpg', dpi=1200)

plt.show()

# %%
merged_df[(merged_df["AAF_Total"]>=0.1) & (merged_df["integration_variant"]!="PhaCinBetalike")][["tseq","integration_variant","Genes","AAF_Total","Gene","Consequence"]].sort_values('AAF_Total',ascending=False)

# %%
merged_df['Gene']=merged_df['Gene'].fillna("")
merged_df[merged_df['Gene'].str.contains('SLC29A1')]

# %%
#merged_df['Genes']=merged_df['Genes'].fillna("")
merged_df['AAF_Total']=merged_df[samples_minusEU].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
merged_df[merged_df['Genes'].str.contains('SLC29A1')][["tseq","start","end","AAF_Total","Genes","CHROMPOS","Consequence","Gene"]]

# %%
result_df['Genes']=result_df['Genes'].fillna("")
result_df[result_df["Genes"].str.contains('SLC29A1')]

# %%
merged_df.loc[370,'INFO_y']

# %%
sns.set_theme(style="ticks")
aaf2 = pd.DataFrame()
aaf2['Alternative Allele Frequency']=merged_df[merged_df["AAF_NSW"]>0]['AAF_NSW']
aaf2['State']='NSW'
aaf2['integration_variant']=merged_df[merged_df["AAF_NSW"]>0]['integration_variant']
aaf = pd.DataFrame()
aaf['Alternative Allele Frequency']=merged_df[merged_df["AAF_QLD"]>0]['AAF_QLD']
aaf['State']='QLD'
aaf['integration_variant']=merged_df[merged_df["AAF_QLD"]>0]['integration_variant']
aaf=pd.concat([aaf,aaf2],ignore_index=True)
aaf3 = pd.DataFrame()
aaf3['Alternative Allele Frequency']=merged_df[merged_df["AAF_VIC"]>0]['AAF_VIC']
aaf3['State']='VIC'
aaf3['integration_variant']=merged_df[merged_df["AAF_VIC"]>0]['integration_variant']
aaf=pd.concat([aaf,aaf3],ignore_index=True)

plt.figure(figsize=(10, 6))
g = sns.boxplot(x="integration_variant", y="Alternative Allele Frequency",
            data=aaf.sort_values("integration_variant"), palette=["mistyrose","orange","lightblue"],hue='State',fliersize=2,linewidth=2)
locs, labels = plt.xticks() 
plt.xticks(locs,["enKoRV \n(n=2,616)",r'phaCin-$\beta$'+' \n(n=1,446)',r'phaCin-$\beta$-like'+' \n(n=359)'])

plt.subplots_adjust(bottom=0.3)  # Adjust the bottom to make room for the legend
plt.xlabel('')
# Place the legend below the plot
plt.legend(loc='center', bbox_to_anchor=(0.45, -0.3), ncol=3, borderaxespad=0., frameon=False)
sns.despine(offset=10, trim=True)
plt.savefig('AAF_wild_IS.jpg', dpi=1200)
plt.show()

# %%
merged_df[(merged_df["AAF_VIC"]>0.4) & (merged_df["integration_variant"]=="KoRV")]

# %%
vcf_df[vcf_df['Gene']=='SLC12A7']


