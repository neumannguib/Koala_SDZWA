# %%
import glob
import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.formula.api import logit
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score
from scipy.stats import mode
import os
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')
sns.set_theme(style="ticks", palette="colorblind")

# %%
#get sample info 
path = 'path to Koala_Genome_Survey/'
sample_info = pd.read_excel("samples_selected.xlsx")
sample_info.groupby('State').count()['Sample ID']

# %%
df_IS = pd.DataFrame()
for sample in sample_info['AWS File Name']:
    try:
        temp = pd.read_csv(path+str(sample)+"_integration_sites.tsv",sep='\t')
        temp['Sample'] = sample
        df_IS = pd.concat([df_IS,temp])
    except:
        print(sample)
df_IS

# %%
df_IS.describe()

# %%
len(df_IS['Sample'].unique()) #number of koala samples

# %%
samples = df_IS['Sample'].unique()
df_samples = pd.DataFrame(samples,columns=['Sample'])
df_samples = pd.merge(df_samples, sample_info, how='left', left_on='Sample', right_on='AWS File Name')
df_samples

# %%
df_IS.groupby('integration_variant').count()['tstart'] #total number of IS per type for all samples

# %%
df_IS['windown_len']=abs(df_IS['tend']-df_IS['tstart'])
df_IS.describe()

# %%
variants = df_IS.groupby(by=['Sample','integration_variant']).count()['tstart'].reset_index()
variants.groupby('integration_variant')['tstart'].max()

# %%
variants.groupby('integration_variant')['tstart'].sum()/len(df_IS['Sample'].unique()) 

# %%
len(df_IS)

# %%
# Sort the DataFrame by POSITIONS + TSD
df_IS=df_IS[~df_IS['POS'].isna()].reset_index(drop=True)#use only IS with POS defined
TSD=10
df_IS['start']=df_IS['POS'].apply(lambda x: int(x.split('-')[0])-TSD if '-' in str(x) else int(x)-TSD)
df_IS['end']=df_IS['POS'].apply(lambda x: int(x.split('-')[1])+TSD if '-' in str(x) else int(x)+TSD)
df_IS = df_IS.sort_values(by=['integration_variant','tseq','start']).reset_index(drop=True)

# Initialize a column to keep track of merged groups
df_IS['group'] = range(1, len(df_IS) + 1)  # Each row starts in its own group

# Function to find overlapping groups and update the 'group' column based on tstart and tend of clusters 
def update_groups(row):
    current_group = row['group']
    overlapping_groups = df_IS[(df_IS['tseq']==row['tseq']) & (df_IS['integration_variant']==row['integration_variant']) & (df_IS['end'] >= row['start']) & 
    (df_IS['start'] <= row['end']) & 
        (df_IS['group'] != current_group)]['group'].tolist()
    if overlapping_groups:
        min_group = min([current_group] + overlapping_groups)
        df_IS.loc[df_IS['group'].isin([current_group] + overlapping_groups), 'group'] = min_group

# Iterate through rows and update the 'group' column
df_IS.apply(update_groups, axis=1)
print(len(df_IS))

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
'end': get_mode,'POS': [('POS', get_mode), ('POS_mode_count', get_mode_count)],'Sample':'nunique','reads_evidence':'mean',"group":"first"}).reset_index(drop=True)


# %%
result_df.columns = result_df.columns.get_level_values(0)
result_df

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
result_df['Freq']=result_df['Sample']/len(df_IS['Sample'].unique()) 
result_df.sort_values('Freq',ascending=False)

# %%
result_df.groupby(by='integration_variant').count()['start']

# %%
result_df[(result_df['Freq']>0.4) & (result_df["integration_variant"]=="KoRV")]

# %%
# Draw a nested boxplot to show bills by day and time
sns.boxplot(x="integration_variant", y="Freq", data=result_df)
sns.despine(offset=10, trim=True)

# %%
#masked info - > check which integrations
masked = pd.read_csv("GCA_030178435.1_masked_KoRV_PhER_PhaCinBeta/GCA_030178435.1_ASM3017843v1_genomic.fna.out",
sep='\s+',skiprows=3,names=['score','div.','del.','ins.','sequence','tstart','tend','left','.','repeat','class/family','begin','end','(left)','ID','extra'])
masked=masked[masked['repeat'].str.contains("PhER|PhaCinBeta|AB721500.1|AF151794.2|KC779547.1")].reset_index(drop=True)
masked['len']=masked['tend']-masked['tstart']
masked['repeat']=masked['repeat'].replace({"AB721500.1":"KoRV","AF151794.2":"KoRV","KC779547.1":"KoRV"})
masked['repeat']=masked['repeat'].apply(lambda x: "PhaCinBetalike" if "Betalike" in x else "PhaCinBeta" if "Beta" in x else x)
masked=masked[masked['len']>10].reset_index(drop=True)
is_masked=[]
for i in df_IS.index:
    if "-" in str(df_IS.loc[i,'POS']):
        pos1=int(df_IS.loc[i,'POS'].split("-")[0]);pos2=int(df_IS.loc[i,'POS'].split("-")[1])
        if len(masked[(masked['sequence']==df_IS.loc[i,"tseq"]) & ((abs(masked['tstart']-pos1)<=TSD)
        & (abs(masked['tend']-pos2)<=TSD))])>0:
            is_masked.append(True)
            continue
    is_masked.append(False)
df_IS['masked']=is_masked
df_IS['masked'].sum()

# %%
df_IS[df_IS['masked']].groupby('integration_variant')['group'].describe()

# %%
df_IS[~df_IS['masked']].groupby('integration_variant')['group'].describe()

# %%
masked.groupby(['repeat','class/family']).count()['ID']

# %%
masked.groupby('repeat')['len'].describe()

# %%
df_IS.to_csv(path+"all_IS_Genome_Survey.tab",sep="\t",index=False)

# %%
result_df.to_csv(path+"IS_frequencies_Genome_Survey.tab",sep="\t",index=False)
result_df.describe()

# %%
#get results from zygosity (Step7-2)
result_df2= pd.read_csv(path+"IS_frequencies_Genome_Survey.tab",sep="\t")
result_df2

# %%
sample_info=sample_info.astype({'AWS File Name':str})
samples = list(sample_info['AWS File Name'].unique()) #taking the duplicates consensus 
result_df2.loc[:,samples].fillna(0)
result_df2['#Missing']= result_df2[samples].apply(lambda row: (row < 0).sum(), axis=1)
result_df2['AAF']= result_df2[samples].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
result_df2['Heterozygosity']=result_df2[samples].apply(lambda row: (row == 1 ).sum()/(row >=0 ).sum(), axis=1)
result_df2['Homozygous']=result_df2[samples].apply(lambda row: (row == 2 ).sum()/(row >=0 ).sum(), axis=1)
result_df2.describe()

# %%
result_df2.groupby("integration_variant").count()['group']

# %%
result_df2['Sample']=result_df2[samples].apply(lambda row: (row != 0 ).sum(), axis=1)
result_df2['Missing']=result_df2['#Missing']/75
result_df2[(result_df2['AAF']>0) & (result_df2['Missing']<=0.1)].groupby("integration_variant").count()['group']


# %%
result_df3=result_df2[(result_df2['AAF']>0) & (result_df2['Missing']<=0.1) ]
result_df3.groupby("integration_variant").count()['group']

# %%
result_df2[(result_df2['AAF']<1) & (result_df2['AAF']>0) & (result_df2['Missing']<=0.1)].groupby("integration_variant").count()['group']

# %%
result_df3[result_df3["integration_variant"]=="PhaCinBetalike"].describe()

# %%
plt.figure(figsize=(10, 6))
g = sns.boxplot(x="integration_variant", y="Homozygous",
            data=result_df3)
g.set_title('Homozygous Insertion')
locs, labels = plt.xticks()
plt.xticks(locs,["KoRV \n(n=1,616)","PhaCinBeta \n(n=1,445)","PhaCinBetalike \n(n=288)"])
sns.despine(offset=10, trim=True)
plt.show()

# %%
count=[]; reference_length=3234982288 #wilpena
for sample in sample_info['Sample ID_dup']:
    with open('path to alignments/'+str(sample)+'_sorted.bam_stats.txt', 'r') as file:
        for line in file:
            if "bases mapped:" in line:  # Find total bases mapped
                total_mapped_bases = int(line.split(':')[1].split('#')[0].strip())
                break
        count.append(total_mapped_bases / reference_length)
sample_info['Coverage']=count
sample_info.describe()
# %%
result_df3.to_excel("IS_all_wild.xlsx",index=False)


