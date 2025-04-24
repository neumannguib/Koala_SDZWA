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
path = 'path to processed data/'
sample_info = pd.read_excel('metadata.xlsx')

# %%
count=[]; reference_length=3234982288 #wilpena
for sample in sample_info['Sample ID_dup']:
    with open(path+str(sample)+'_sorted.bam_stats.txt', 'r') as file:
        for line in file:
            if "bases mapped:" in line:  # Find total bases mapped
                total_mapped_bases = int(line.split(':')[1].split('#')[0].strip())
                break
        count.append(total_mapped_bases / reference_length)
sample_info['Coverage']=count
sample_info.describe()

# %%
df_IS = pd.DataFrame()
for sample in sample_info['Sample ID_dup']:
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
df_IS.groupby('integration_variant').count()['tstart'] #total number of IS per type for all samples

# %%
variants = df_IS.groupby(by=['Sample','integration_variant']).count()['tstart'].reset_index()
variants.groupby('integration_variant')['tstart'].max()

# %%
variants.groupby('integration_variant')['tstart'].sum()/len(df_IS['Sample'].unique()) 

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
# Draw a nested boxplot to show bills by day and time
sns.boxplot(x="integration_variant", y="Freq", data=result_df)
sns.despine(offset=10, trim=True)

# %%
#masked info - > check which integrations
TSD=10
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
result_df.to_csv(path+"IS_frequencies.tab",sep="\t",index=False)
result_df.describe()

# %%
df_IS.to_csv(path+"all_IS.tab",sep="\t",index=False)

# %%
#get results with zygosity
result_df2= pd.read_csv(path+"IS_frequenciesv2.tab",sep="\t")
result_df2

# %%
#get consensus for replicates
replicates=['153','169','308']
result_df2['153']=[result_df2.loc[i,'153_liver'] if result_df2.loc[i,'153_liver']==result_df2.loc[i,'153_lung'] else -1 for i in result_df2.index]
result_df2['169']=[result_df2.loc[i,'169_nose'] if result_df2.loc[i,'169_nose']==result_df2.loc[i,'169_lung'] else -1 for i in result_df2.index]
result_df2['308']=[result_df2.loc[i,'308_liver'] if result_df2.loc[i,'308_liver']==result_df2.loc[i,'308_spleen'] else -1 for i in result_df2.index]

# %%
print('153 '+str((result_df2['153_liver']==result_df2['153_lung']).sum()/len(result_df2)))
print('169 '+str((result_df2['169_nose']==result_df2['169_lung']).sum()/len(result_df2)))
print('308 '+str((result_df2['308_liver']==result_df2['308_spleen']).sum()/len(result_df2)))

# %%
print('153 '+str(len(result_df2[(result_df2['153_lung']>=0)&(result_df2['153_liver']>=0)&(result_df2['153_liver']==result_df2['153_lung'])])/len(result_df2[(result_df2['153_lung']>=0)&(result_df2['153_liver']>=0)])))
print('169 '+str(len(result_df2[(result_df2['169_lung']>=0)&(result_df2['169_nose']>=0)&(result_df2['169_nose']==result_df2['169_lung'])])/len(result_df2[(result_df2['169_lung']>=0)&(result_df2['169_nose']>=0)])))
print('308 '+str(len(result_df2[(result_df2['308_spleen']>=0)&(result_df2['308_liver']>=0)&(result_df2['308_liver']==result_df2['308_spleen'])])/len(result_df2[(result_df2['308_spleen']>=0)&(result_df2['308_liver']>=0)])))

# %%
print('153 '+str(len(result_df2[(result_df2['integration_variant']=="KoRV")&(result_df2['153_liver']==result_df2['153_lung'])])/len(result_df2[(result_df2['integration_variant']=="KoRV")])))
print('169 '+str(len(result_df2[(result_df2['integration_variant']=="KoRV")&(result_df2['169_nose']==result_df2['169_lung'])])/len(result_df2[(result_df2['integration_variant']=="KoRV")])))
print('308 '+str(len(result_df2[(result_df2['integration_variant']=="KoRV")&(result_df2['308_liver']==result_df2['308_spleen'])])/len(result_df2[(result_df2['integration_variant']=="KoRV")])))

# %%
print(len(result_df2[(result_df2['153_lung']>=0)&(result_df2['153_liver']>=0)]))
print(len(result_df2[(result_df2['169_lung']>=0)&(result_df2['169_nose']>=0)]))
print(len(result_df2[(result_df2['308_spleen']>=0)&(result_df2['308_liver']>=0)]))

# %%
result_df2[result_df2['153_liver']!=result_df2['153_lung']].groupby('integration_variant').count()['group']

# %%
result_df2[(result_df2['153_liver']!=result_df2['153_lung']) & (result_df2['integration_variant']=='KoRV')][['tseq','start','end','POS','153_liver','153_lung']]

# %%
sample_info=sample_info.astype({'Sample ID':str})
samples = list(sample_info['Sample ID'].unique()) #taking the duplicates consensus 
result_df2.loc[:,samples].fillna(0)
result_df2['#Missing']= result_df2[samples].apply(lambda row: (row < 0).sum(), axis=1)
#result_df2['#Somatic']= result_df2[samples].apply(lambda row: (row == -2).sum(), axis=1)
result_df2['AAF']= result_df2[samples].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
result_df2['Freq']=result_df2[samples].apply(lambda row: (row != 0 ).sum()/len(samples), axis=1)
result_df2['Heterozygosity']=result_df2[samples].apply(lambda row: (row == 1 ).sum()/(row >=0 ).sum(), axis=1)
result_df2['Homozygous']=result_df2[samples].apply(lambda row: (row == 2 ).sum()/(row >=0 ).sum(), axis=1)
result_df2.describe()

# %%
result_df2.groupby("integration_variant").count()['group']

# %%
result_df2['Sample']=result_df2[samples].apply(lambda row: (row != 0 ).sum(), axis=1)
result_df2['Missing']=result_df2['#Missing']/91
#result_df2['Somatic']=result_df2['#Somatic']/91
result_df2[(result_df2['AAF']>0) & (result_df2['Missing']<=0.1)].groupby("integration_variant").count()['group']


# %%
pheno = pd.read_excel("Koala health outcome data.xlsx",skiprows=1)
pheno = pheno.rename({"SB":"IndividualId","Birth_location_c":"FamilyId","Sex_c":'Gender',"Sire":"Father","Dam":"Mother","Dead":"Deceased"},axis=1)
pheno = pheno[pheno["WGS_batch"]!="Fail"].reset_index(drop=True)
pheno["Deceased"]=pheno["Deceased"].replace({"Dead":"Y","Alive":"."})
pheno["Gender"]=pheno["Gender"].replace({"Male":"M","Female":"F"})
pheno["Father"]=pheno["Father"].replace({9999:"."})
pheno['Mother']=pheno['Mother'].replace({9999:"."})
pheno["COD3"]=pheno["COD3"].replace({9999:""})
pheno.index=pheno['IndividualId']

# %%
pheno

# %%
#count how many trios
c=0;fathers=[];mothers=[]
f1=[];f2=[];f3=[];f4=[]
for index, row in pheno.sort_values('Accession').iterrows():
    sample, father, mother = row['Accession'], row['Father'], row['Mother']
    father_genotyped = father in pheno['IndividualId']
    mother_genotyped = mother in pheno['IndividualId']
    if father_genotyped:
        father=pheno.loc[father,'Accession']
        fathers.append(father)
    if mother_genotyped:
        mother=pheno.loc[mother,'Accession']
        mothers.append(mother)
    if father_genotyped and mother_genotyped:
        c+=1        
        if father in f3 or mother in f3:
            f4.append(sample)
        elif father in f2 or mother in f2:
            f3.append(sample)
        elif father in f1 or mother in f1:
            f2.append(sample)
        else:
            f1.append(sample)
        print(sample, father, mother)

# %%
def check_mendelian_error_with_single_parent(child_gt, father_gt, mother_gt):
    """
    Check for Mendelian errors given genotypes of a child and their parents.
    Now also handles cases with only one parent's genotype available.
    Genotypes are encoded as 0, 1, 2, with -1 and -2 representing missing data.
    """
    # Define missing data handling
    if child_gt <0:
        return np.nan  # Inconclusive due to missing data in child

    # Scenario with both parents
    if father_gt >=0 and mother_gt >=0:
        if (father_gt == 0 and mother_gt == 2 and child_gt !=1) or (father_gt == 2 and mother_gt == 0 and child_gt != 1):
            return True  # Mendelian error
        if (father_gt in [0, 2] and mother_gt in [0, 2] and father_gt == mother_gt and child_gt != father_gt):
            return True  # Homozygous parents but child's genotype doesn't match
        return False  # No error detected

    # Scenario with only one parent
    parent_gt = father_gt if father_gt >=0 else mother_gt
    if parent_gt in [0, 2]:
        if (parent_gt == 0 and child_gt == 2) or (parent_gt == 2 and child_gt == 0):
            return True  # Mendelian error for opposing homozygous
    return False  # No error detected if one parent is heterozygous or data is missing


errors = []
samples_count={}
for index, row in pheno.iterrows():
    sample, father, mother = row['Accession'], row['Father'], row['Mother']
    father_genotyped = father in pheno['IndividualId']
    mother_genotyped = mother in pheno['IndividualId']

    if father_genotyped or mother_genotyped:
        if father_genotyped:
            father=str(pheno.loc[father,'Accession'])
        if mother_genotyped:
            mother=str(pheno.loc[mother,'Accession'])
        samples_count[sample]=0
        for variant in result_df2.index:
            
            child_gt = result_df2.at[variant, str(sample)]
            if child_gt in [-1, -2]:
                continue
            father_gt = result_df2.at[variant, father] if father_genotyped else -1  # Use -1 for missing data
            mother_gt = result_df2.at[variant, mother] if mother_genotyped else -1

            if check_mendelian_error_with_single_parent(child_gt, father_gt, mother_gt):
                samples_count[sample]+=1
                integration_variant=result_df2.at[variant, 'integration_variant']
                group=result_df2.at[variant, 'group']
                count_is=result_df2.at[variant, 'Sample'] 
                variant=result_df2.at[variant, 'tseq']+":"+result_df2.at[variant, 'POS']
                errors.append((sample, variant, group, integration_variant, child_gt, father, father_gt, mother, mother_gt,count_is))

# Convert errors to a DataFrame
errors_df = pd.DataFrame(errors, columns=['Sample', 'Variant',"group",'Type', 'Child_GT', 'Father_ID','Father_GT', 'Mother_ID', 'Mother_GT','Count'])
errors_df

# %%
errors_df['Child_GT'].unique()

# %%
temp = errors_df[(errors_df['Child_GT']==1) & (errors_df['Father_GT']==0) & (errors_df['Mother_GT']==0)].reset_index(drop=True)
covs=[];cov_father=[];cov_mother=[];f=[]
temp=temp.astype({'Sample':str,'group':int})
df_IS=df_IS.astype({'Sample':str,'group':int})

for i in temp.index:
    covs.append(int(df_IS[(df_IS['Sample']==temp.loc[i,'Sample']) & (df_IS['group']==temp.loc[i,'group'])]['reads_evidence'].iloc[0]))
    cov_father.append(df_IS[(df_IS['Sample']==temp.loc[i,'Father_ID'])]['reads_evidence'].mean())
    cov_mother.append(df_IS[(df_IS['Sample']==temp.loc[i,'Mother_ID'])]['reads_evidence'].mean())
    sample=int(temp.loc[i,'Sample'])
    f.append( fathers.count(sample) + mothers.count(sample))
temp['Coverage_child']=covs
temp['Coverage_father']=cov_father
temp['Coverage_mother']=cov_mother
temp['fx']=f
temp

# %%
result_df2.index=result_df2['group']
#validate one by one on the f1 and f2
pheno.index=pheno['Accession']
def validate_offspring(sample,group):
    if sample in fathers:
        father=pheno.loc[sample,'IndividualId']
        temp_pheno = pheno[pheno['Father']==father]
        c=0
        for code in temp_pheno['Accession']:
            if result_df2.loc[group,str(code)]>0:
                if code in fathers or code in mothers:
                    c+=1 + validate_offspring(code,group)
                else:
                    c+=1
        return c
    mother=pheno.loc[sample,'IndividualId']
    temp_pheno = pheno[pheno['Mother']==mother]
    c=0
    for code in temp_pheno['Accession']:
        if result_df2.loc[group,str(code)]>0:
            if code in fathers or code in mothers:
                c+=1 + validate_offspring(code,group)
            else:
                c+=1
    return c

validate=[]
for i,sample in enumerate(temp['Sample']):
    group=temp.loc[i,'group']
    sample = int(sample)
    if temp.loc[i,'Count']>1:
        if sample in fathers or sample in mothers:
            if temp.loc[i,'Count'] == (1+validate_offspring(sample,group)):
                validate.append(True)
            else:
                validate.append(False)
        else:
            validate.append(False)
    else:
        validate.append(True)
temp['Validation_offspring']=validate
temp[temp['Validation_offspring']]

# %%
temp[(temp['Validation_offspring']) & (temp['Coverage_child']>=40)].groupby(['Sample','Type'])['Variant'].count()

# %%
temp[(temp['Validation_offspring']) & (temp['Coverage_child']>=40)]

# %%
temp=pd.merge(temp.astype({'Sample':int}),pheno.reset_index(drop=True)[['Accession','Gender','FamilyId','DOB', 'DOD', 'Age1','COD_dx','Sire_age','Dam_age']],how='left',left_on='Sample',right_on='Accession')
temp[(temp['Validation_offspring']) & (temp['Coverage_child']>=40)].drop_duplicates('Sample',keep='first')[['Sample','Gender','FamilyId','DOB', 'DOD', 'Age1','COD_dx','Sire_age','Dam_age']]

# %%
temp=pd.merge(temp.astype({'Sample':int}),sample_info.astype({'Sample ID':int})[['Sample ID','Tissue Source','Date of collection','Geographical Location']],how='left',left_on='Sample',right_on='Sample ID')
temp[(temp['Validation_offspring']) & (temp['Coverage_child']>=10)].drop_duplicates('Sample',keep='first')[['Sample','Gender','FamilyId','DOB', 'DOD', 'Age1','COD_dx','Sire_age','Dam_age','Tissue Source','Date of collection','Geographical Location']]

# %%
total_sequence_length = 2*(3234982288) #wilpena
15 / (total_sequence_length * 46)

# %%
total_sequence_length = 2*(3234982288) #wilpena
1 / (total_sequence_length * 46)

# %%
errors_df.to_csv(path+"mendelian_errors.tab",sep="\t",index=False)

# %%
temp.to_csv(path+"denovoERVs.tab",sep="\t",index=False)

# %%
from scipy.stats import binom
cum_pop=0
for reads in range(12,29):
    p=binom.pmf(reads, 40, 0.5)
    print(p)
    cum_pop+= p
print(1-cum_pop)

# %%
from scipy.stats import binom
cum_pop=0
for reads in range(3,8):
    p=binom.pmf(reads, 10, 0.5)
    print(p)
    cum_pop+= p
print(1-cum_pop)

# %%
cum_pop

# %%
cum_pop=0
for reads in range(0,6):
    p=binom.pmf(reads, 16, 0.45)
    print(p)
    cum_pop+= p
cum_pop

# %%
df_erro_perIS = errors_df.groupby(by="group").count()["Sample"].to_frame(name="Mendel_error").reset_index()
result_df2 = pd.merge(result_df2.reset_index(drop=True), df_erro_perIS, how="left", on="group" )
result_df2["Mendel_error"] = result_df2["Mendel_error"].fillna(0)
result_df2["Mendel_error"] = result_df2["Mendel_error"] / 63
result_df2["Mendel_error"].describe()

# %%
result_df3=result_df2[(result_df2['AAF']>0) & (result_df2['Missing']<=0.1) & (result_df2["Mendel_error"]<=0.05)]
result_df3.groupby("integration_variant").count()['group']

# %%
result_df2[(result_df2['AAF']<1) & (result_df2['AAF']>0) & (result_df2['Missing']<=0.1) & (result_df2["Mendel_error"]<=0.05)].groupby("integration_variant").count()['group']

# %%
result_df3[result_df3["integration_variant"]=="PhaCinBetalike"].describe()

# %%
pheno['KoRV']=[len(result_df3[(result_df3['integration_variant']=='KoRV') & (result_df3[str(sample)]>0)]) for sample in pheno['Accession']]
pheno['phaCinBeta']=[len(result_df3[(result_df3['integration_variant']=='PhaCinBeta') & (result_df3[str(sample)]>0)]) for sample in pheno['Accession']]
pheno['phaCinBetalike']=[len(result_df3[(result_df3['integration_variant']=='PhaCinBetalike') & (result_df3[str(sample)]>0)]) for sample in pheno['Accession']]
pheno['Total IS']=pheno['KoRV']+pheno['phaCinBeta']+pheno['phaCinBetalike']
pheno[['Total IS','KoRV','phaCinBeta','phaCinBetalike']].describe()

# %%
#calculate mendelian error per animal after filtering
errors_dfv2=errors_df[errors_df['Type']!='PhER'].reset_index(drop=True)
remove=[]
for i,group in enumerate(errors_dfv2['group']):
    if group not in list(result_df3['group']):
        remove.append(i)
errors_dfv2=errors_dfv2.drop(remove,axis=0).reset_index(drop=True)
mendelian_error=pd.DataFrame()
mendelian_error['Sample']=samples_count.keys()
pheno.index=pheno['Accession']
mendelian_error['n_variants']=list(pheno.loc[list(mendelian_error['Sample']),"Total IS"])
mendelian_error=mendelian_error.sort_values('Sample')
mendelian_error['n_errors']=[len(errors_dfv2[errors_dfv2['Sample']==sample]) for sample in mendelian_error['Sample']]
mendelian_error['Error rate (%)']=mendelian_error['n_errors']/mendelian_error['n_variants']*100

mendelian_error.index=mendelian_error['Sample']
pheno['Mendelian error(%)']=[mendelian_error.loc[sample,'Error rate (%)'] if sample in list(mendelian_error['Sample']) else 0 for sample in pheno['Accession']]
pheno.describe()#but only 63 had parents, so mean from 63 should be used

mendelian_error.describe()

# %%
mendelian_error.to_csv(path+"mendelian_errors_per_sample.tab",sep="\t",index=False)

# %%
mendelian_error[mendelian_error['Error rate (%)']>1]

# %%
result_df3.groupby("integration_variant").count()['group']

# %%
len(result_df3)

# %%
plt.figure(figsize=(10, 6))
g = sns.boxplot(x="integration_variant", y="Homozygous",
            data=result_df3)
g.set_title('Homozygous Insertion')
locs, labels = plt.xticks()
plt.xticks(locs,["KoRV \n(n=1,232)","PhaCinBeta \n(n=485)","PhaCinBetalike \n(n=358)"])
sns.despine(offset=10, trim=True)
plt.savefig(path+'homozygous_ISv4.jpg', dpi=600)
plt.show()

# %%
count=[]; reference_length=3234982288 #wilpena
for sample in sample_info['Sample ID_dup']:
    with open(path+str(sample)+'_sorted.bam_stats.txt', 'r') as file:
        for line in file:
            if "bases mapped:" in line:  # Find total bases mapped
                total_mapped_bases = int(line.split(':')[1].split('#')[0].strip())
                break
        count.append(total_mapped_bases / reference_length)
sample_info['Coverage']=count
sample_info.describe()

# %%
c=[]
for sample in pheno["Accession"]:
    if len(sample_info[sample_info["Sample ID"]==str(sample)])>1:
        c.append(sample_info[sample_info["Sample ID"]==str(sample)]['Coverage'].mean())
    else:
        c.append(list(sample_info[sample_info["Sample ID"]==str(sample)]['Coverage'])[-1])
pheno['Coverage']=c
pheno[pheno['Coverage']<30]

# %%
pheno["Birth_year"]=pheno["DOB"].apply(lambda x: x.year)
g = sns.lmplot(
    data=pheno,
    x="Birth_year", y="KoRV", hue="Deceased",
    height=5)

# %%
g = sns.lmplot(
    data=pheno,
    x="Birth_year", y="phaCinBeta", hue="Deceased",
    height=5)

# %%
g = sns.lmplot(
    data=pheno,
    x="Birth_year", y="phaCinBetalike", hue="Deceased",
    height=5
)

# %%
samples_alive = [str(sample) for sample in pheno[pheno['Dead_c']!="Dead"]['Accession']]
samples_dead = [str(sample) for sample in pheno[pheno['Dead_c']=="Dead"]['Accession']]
result_df3['AAF-alive']= result_df3[samples_alive].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
result_df3['AAF-dead']= result_df3[samples_dead].apply(lambda row: ((row == 1 ).sum()+((row == 2).sum()*2))/((row >= 0).sum()*2), axis=1)
result_df3['AAF-diff'] = result_df3['AAF-alive'] - result_df3['AAF-dead']
result_df3.describe()

# %%
result_df3[(result_df3['AAF-alive']<1) & (result_df3['integration_variant']=='PhaCinBetalike')]

# %%
sns.set_theme(style="ticks")
aaf2 = pd.DataFrame()
aaf2['Alternative Allele Frequency']=result_df3['AAF-alive']
aaf2['Status']='alive (n=33)'
aaf2['integration_variant']=result_df3['integration_variant']
aaf = pd.DataFrame()
aaf['Alternative Allele Frequency']=result_df3['AAF-dead']
aaf['Status']='deceased (n=58)'
aaf['integration_variant']=result_df3['integration_variant']
aaf=pd.concat([aaf,aaf2],ignore_index=True)
plt.figure(figsize=(10, 6))
g = sns.boxplot(x="integration_variant", y="Alternative Allele Frequency",
            data=aaf,palette=['darkgray','sandybrown'], hue='Status',fliersize=2,linewidth=2)
locs, labels = plt.xticks() 
plt.xticks(locs,["KoRV \n(n=1,232)",r'phaCin-$\beta$'+' \n(n=485)',r'phaCin-$\beta$-like'+' \n(n=358)'])

plt.subplots_adjust(bottom=0.3)  # Adjust the bottom to make room for the legend
plt.xlabel('')
# Place the legend below the plot
plt.legend(loc='center', bbox_to_anchor=(0.45, -0.3), ncol=2, borderaxespad=0., frameon=False)
sns.despine(offset=10, trim=True)
plt.savefig(path+'AAF_alive_deadv5.jpg', dpi=1200)
plt.show()

# %%
sns.set_theme(style="ticks")
aaf2 = pd.DataFrame()
aaf2['Alternative Allele Frequency']=result_df3['AAF-alive']
aaf2['Status']='alive (n=33)'
aaf2['integration_variant']=result_df3['integration_variant']
aaf = pd.DataFrame()
aaf['Alternative Allele Frequency']=result_df3['AAF-dead']
aaf['Status']='deceased (n=58)'
aaf['integration_variant']=result_df3['integration_variant']
aaf=pd.concat([aaf,aaf2],ignore_index=True)
plt.figure(figsize=(10, 6))
g = sns.boxplot(x="integration_variant", y="Alternative Allele Frequency",
            data=aaf,palette=['darkgray','forestgreen'], hue='Status',fliersize=2,linewidth=2)
locs, labels = plt.xticks() 
plt.xticks(locs,["enKoRV \n(n=1,232)",r'phaCin-$\beta$'+' \n(n=485)',r'phaCin-$\beta$-like'+' \n(n=358)'])

plt.subplots_adjust(bottom=0.3)  # Adjust the bottom to make room for the legend
plt.xlabel('')
# Place the legend below the plot
plt.legend(loc='center', bbox_to_anchor=(0.45, -0.3), ncol=2, borderaxespad=0., frameon=False)
sns.despine(offset=10, trim=True)
plt.show()

# %%
path

# %%
fig = plt.figure(figsize=(10, 4))
ax = sns.boxplot(x="integration_variant", y="Alternative Allele Frequency",
            data=aaf,hue='Status',palette="Set2",fliersize=5,linewidth=1,linecolor='white')
locs, labels = plt.xticks() 
# Set axes spines to white
ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')
ax.spines['bottom'].set_color('white')
ax.spines['left'].set_color('white')

# Set axes ticks to white
ax.xaxis.label.set_color('white')
ax.yaxis.label.set_color('white')
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')

# Set figure background to transparent and axes background to black (to see white lines)
fig.patch.set_alpha(0.0)
plt.xticks(locs,["KoRV \n(n=1,232)",r'phaCin-$\beta$'+' \n(n=485)',r'phaCin-$\beta$-like'+' \n(n=358)'])
ax.tick_params(axis='both', which='major', labelsize=20)
ax.tick_params(axis='both', which='minor', labelsize=20)
ax.set_xlabel('', fontsize=12)
ax.set_ylabel('Insertion Allele Frequency', fontsize=15)
sns.despine(offset=10, trim=True)
legend = ax.legend(loc='upper left', bbox_to_anchor=(1.1, 1), borderaxespad=0.,frameon=False,fontsize=20)
# Set legend text color to white
for text in legend.get_texts():
    text.set_color('white')
plt.show()

# %%
plt.figure(figsize=(10, 6))
g = sns.boxplot(x="integration_variant", y="Heterozygosity",
            data=result_df3)
g.set_title('Heterozygosity')
locs, labels = plt.xticks()
plt.xticks(locs,["KoRV \n(n=1,232)",r'phaCin-$\beta$'+' \n(n=485)',r'phaCin-$\beta$-like'+' \n(n=358)'])
sns.despine(offset=10, trim=True)
plt.show()

# %%
plt.figure(figsize=(10, 6))
g = sns.boxplot(x="integration_variant", y="Freq",
            data=result_df2[result_df2['integration_variant']!='PhER'])
g.set_title('Frequency (simple presence x absence)')
locs, labels = plt.xticks()
plt.xticks(locs,["KoRV \n(n=1,232)",r'phaCin-$\beta$'+' \n(n=485)',r'phaCin-$\beta$-like'+' \n(n=358)'])
sns.despine(offset=10, trim=True)
plt.show()

# %%
plt.figure(figsize=(10, 6))
g = sns.boxplot(x="integration_variant", y="Freq",
            data=result_df2)
g.set_title('Frequency')
locs, labels = plt.xticks()
#plt.xticks(locs,["KoRV \n(n=1,583)","PhaCinBeta \n(n=628)","PhaCinBetalike \n(n=955)"])
sns.despine(offset=10, trim=True)
plt.show()

# %%
plt.figure(figsize=(10, 6))
g = sns.boxplot(x="integration_variant", y="AAF",
            data=result_df3[result_df3['integration_variant']!='PhER'])
g.set_title('AAF')
locs, labels = plt.xticks()
plt.xticks(locs,["KoRV \n(n=1,232)",r'phaCin-$\beta$'+' \n(n=485)',r'phaCin-$\beta$-like'+' \n(n=358)'])
sns.despine(offset=10, trim=True)
plt.show()

# %%
#testing concordance again, but nor with AAF>0 and non Pher variants
print('153 '+str((result_df3['153_liver']==result_df3['153_lung']).sum()/len(result_df3)))
print('169 '+str((result_df3['169_nose']==result_df3['169_lung']).sum()/len(result_df3)))
print('308 '+str((result_df3['308_liver']==result_df3['308_spleen']).sum()/len(result_df3)))

# %%
f, (ax, ax2) = plt.subplots(2, 1, sharex=True)

# plot the same data on both axes
ax.scatter(y=result_df3['reads_evidence'], x=result_df3['Freq'],
            c=result_df3['integration_variant'].astype('category').cat.codes,
            cmap='viridis', alpha=.6, s=10)

g= ax2.scatter(y=result_df3['reads_evidence'], x=result_df3['Freq'],
            c=result_df3['integration_variant'].astype('category').cat.codes,
            cmap='viridis', alpha=.6, s=10)

# zoom-in / limit the view to different portions of the data
ax2.set_ylim(0, 400)  
ax.set_ylim(500, 1500)  
ax2.set_xlim(0, 1.01)  
ax.set_xlim(0, 1.01)  

# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()


d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_ylabel('Average coverage')
ax2.set_xlabel('IS Frequency')
legend_labels = result_df3['integration_variant'].astype('category').unique()
ax2.legend(handles=g.legend_elements()[0], labels=legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left')

plt.show()


# %%
f, (ax, ax2) = plt.subplots(2, 1, sharex=True)

# plot the same data on both axes
ax.scatter(y=result_df3['reads_evidence'], x=result_df3['AAF'],
            c=result_df3['integration_variant'].astype('category').cat.codes,
            cmap='viridis', alpha=.6, s=10)

g= ax2.scatter(y=result_df3['reads_evidence'], x=result_df3['AAF'],
            c=result_df3['integration_variant'].astype('category').cat.codes,
            cmap='viridis', alpha=.6, s=10)

# zoom-in / limit the view to different portions of the data
ax2.set_ylim(0, 400)  
ax.set_ylim(500, 1500)  
ax2.set_xlim(0, 1.01)  
ax.set_xlim(0, 1.01)  

# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()


d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_ylabel('Average coverage')
ax2.set_xlabel('IS AAF')
legend_labels = result_df3['integration_variant'].astype('category').unique()
ax2.legend(handles=g.legend_elements()[0], labels=legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left')

plt.show()


# %%
pheno['KoRV_corrected']=pheno['KoRV']*(pheno['Coverage']/pheno['Coverage'].mean())
pheno['PhaCinBeta_corrected']=pheno['phaCinBeta']*(pheno['Coverage']/pheno['Coverage'].mean())
pheno['PhaCinBetalike_corrected']=pheno['phaCinBetalike']*(pheno['Coverage']/pheno['Coverage'].mean())
df_virus=pd.DataFrame()
pheno["Birth_year"]=[item.year for item in pheno["DOB"]]
df_virus['Corrected Number of IS per Individual']=list(pheno['KoRV_corrected'])+list(pheno['PhaCinBeta_corrected'])+list(pheno['PhaCinBetalike_corrected'])
df_virus['IS type']=['KoRV']*len(pheno)+[r'phaCin-$\beta$']*len(pheno)+[r'phaCin-$\beta$'+'-like']*len(pheno)
df_virus['Birth Year']=list(pheno['Birth_year'])+list(pheno['Birth_year'])+list(pheno['Birth_year'])
g = sns.lmplot(data=df_virus,x="Birth Year", y="Corrected Number of IS per Individual", hue="IS type",height=5,
palette=['indianred',sns.color_palette("husl", 8)[1],'saddlebrown'],legend=False)
plt.subplots_adjust(bottom=0.2)  # Adjust the bottom to make room for the legend
# Place the legend below the plot
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, borderaxespad=0., frameon=False)

# %%
df_virus=pd.DataFrame()
df_virus['Number of IS per individual']=list(pheno['KoRV'])+list(pheno['phaCinBeta'])+list(pheno['phaCinBetalike'])
df_virus['IS type']=['KoRV']*len(pheno)+[r'phaCin-$\beta$']*len(pheno)+[r'phaCin-$\beta$'+'-like']*len(pheno)
df_virus['Birth Year']=list(pheno['Birth_year'])+list(pheno['Birth_year'])+list(pheno['Birth_year'])
g = sns.lmplot(data=df_virus,x="Birth Year", y="Number of IS per individual", hue="IS type",
height=5,palette=['indianred',sns.color_palette("husl", 8)[1],'saddlebrown'],legend=False)
plt.subplots_adjust(bottom=0.2)  # Adjust the bottom to make room for the legend
# Place the legend below the plot
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, borderaxespad=0., frameon=False)

# %%
df_virus=pd.DataFrame()
df_virus['Number of IS per individual']=list(pheno['KoRV'])+list(pheno['phaCinBeta'])
df_virus['IS type']=['KoRV']*len(pheno)+[r'phaCin-$\beta$']*len(pheno)
df_virus['Birth Year']=list(pheno['Birth_year'])+list(pheno['Birth_year'])
g = sns.lmplot(data=df_virus,x="Birth Year", y="Number of IS per individual", hue="IS type",height=5)

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
flanking= 10000
temp['Gene'] = temp["Variant"].apply(
    lambda x: (
        genes.loc[
            (genes["Contig"] == x.split(":")[0]) &
            (genes["start"] <= int(x.split(":")[1])+flanking) &
            (genes["end"] >= int(x.split(":")[1])-flanking)
        ]["Preferred_name"].iloc[0] if (len(genes[(genes["Contig"] == x.split(":")[0]) &
                                                   (genes["start"] <= int(x.split(":")[1])+flanking) &
                                                   (genes["end"] >= int(x.split(":")[1])-flanking)]) > 0)
        else ""
    )
)
temp[temp["Validation_offspring"]]

# %%
#calculate correlation to age at partum

# %%
overlapping_genes = []
flanking= 10000
genes['Preferred_name']=[genes.loc[i,'seed_ortholog'] if '-'==genes.loc[i,'Preferred_name'] else genes.loc[i,'Preferred_name'] for i in genes.index]
# Iterate over each row in contigs_df
for contig_row in result_df3.itertuples(index=False):
    # Filter genes_df to find overlapping genes for the current contig
    overlapping_genes_contig = genes[
        (genes['Contig'] == contig_row.tseq) &
        (genes['start'] <= contig_row.end+flanking) &
        (genes['end'] >= contig_row.start-flanking)
    ]

    # Append the overlapping genes to the list
    overlapping_genes.append(', '.join(overlapping_genes_contig['Preferred_name'].tolist()))
result_df3['Genes']=overlapping_genes
result_df3

# %%
all_genes=[]
file = open(path+"korv_genesv4.txt","w")
for gene in result_df3[result_df3['integration_variant']=='KoRV']['Genes']:
    gene = gene.split(",")
    for item in gene:
        if item.strip() not in all_genes:
            all_genes.append(item.strip())
            file.write(item.strip()+"\n")
file.close()
len(all_genes)

# %%
all_genes=[]
file = open(path+"korv_genes_aaf>0.1v4.txt","w")
for gene in result_df3[(result_df3['integration_variant']=='KoRV') & (result_df3['AAF']>=0.1)]['Genes']:
    gene = gene.split(",")
    for item in gene:
        if item.strip() not in all_genes:
            all_genes.append(item.strip())
            file.write(item.strip()+"\n")
file.close()
len(all_genes)

# %%
all_genes

# %%
all_genes=[]
for gene in result_df3[(result_df3['integration_variant']=='PhaCinBeta')]['Genes']:
    gene = gene.split(",")
    for item in gene:
        if item.strip() not in all_genes:
            all_genes.append(item.strip())
            print(item.strip())

# %%
result_df3[result_df3['Genes'].str.contains('BCL2L1')]

# %%
pedigree = pd.read_csv("Pedigree/koala_san_diego.txt",sep="\t")
pheno['Sequenced']=999900
pedigree = pd.merge(pedigree,pheno[["IndividualId","DOB","Age1","COD3","Sequenced","Birth_year","Accession"]].reset_index(drop=True),how='outer',left_on='Name',right_on="IndividualId")
pedigree['Sequenced']=pedigree['Sequenced'].fillna(999999)
pedigree['COD3']=pedigree['COD3'].replace({53:52,45:53,32:31})
pedigree=pedigree.fillna("").astype({"Family ":str,"Birth_year":str}).astype({"Birth_year":str})
pedigree['Birth']=pedigree[["Family ","Birth_year"]].agg('-'.join, axis=1)
pedigree

# %%
for sample in pedigree[pedigree['COD3']=='']['Name']:
    if sample in list(pheno['Father']):
        col='Father';col2="Sire_COD3"
    elif sample in list(pheno['Mother']):
        col='Mother';col2="Dam_COD3"
    elif sample in list(pheno['PatGSire']):
        col='PatGSire';col2="PGS_COD3"
    elif sample in list(pheno['PatGDam']):
        col='PatGDam';col2="PGD_COD3"
    elif sample in list(pheno['MatGSire']):
        col='MatGSire';col2="MGS_COD3"
    elif sample in list(pheno['MatGDam']):
        col='MatGDam';col2="MGD_COD3"
    else:
        continue
    code = list(pheno[pheno[col]==sample][col2])[-1]
    i = list(pedigree['Name']).index(sample)
    pedigree.loc[i,'COD3']=code
pedigree

# %%
pheno['COD_dx']=pheno['COD_dx'].fillna("");pheno['COD_dx']=pheno['COD_dx'].replace(9999,"")
pheno=pheno.rename({"Accession":"Sample ID"},axis=1)
temp = pheno[pheno['COD_dx'].str.contains("Lymphoma|leukemia|Leukemia|Adenocarcinoma|sarcoma|Osteochondroma|osteochondroma|Sarcoma")]
n = len(temp)
temp

# %%
#Getting freq of IS in koalas which died from cancer
samples_cancer=[str(x) for x in temp['Sample ID']]
result_df3['Freq_cancer_IS'] = [(result_df3.loc[i,samples_cancer]==1).sum()/n for i in result_df3.index]

# %%
#comparing the koalas which died from cancer to other dead koalas 
temp = pheno[(~pheno['COD_dx'].str.contains("Lymphoma|leukemia|Leukemia|Adenocarcinoma|sarcoma|Osteochondroma|osteochondroma|Sarcoma")) & (pheno['Deceased']==1)]
n = len(temp)
samples_notcancer=[str(x) for x in temp['Sample ID']]
result_df3['Freq_death_notcancer_IS'] = [(result_df3.loc[i,samples_notcancer]==1).sum()/n for i in result_df3.index]
result_df3[(result_df3['Freq_cancer_IS']-result_df3['Freq_death_notcancer_IS']>=0.1)
 & (result_df3['integration_variant']=='KoRV')]

# %%
result_df3[(result_df3['Freq_cancer_IS']-result_df3['Freq_death_notcancer_IS']>=0.1)
 & (result_df3['integration_variant']=='PhaCinBeta')]

# %%
result_df3[(result_df3['Freq_cancer_IS']-result_df3['Freq_death_notcancer_IS']>=0.1)][["tseq",	"integration_variant","start",
"end","POS","reads_evidence","Missing",	"AAF-alive","AAF-dead","Genes","Freq_cancer_IS","Freq_death_notcancer_IS"]].to_excel(path+"sup_table1.xlsx",index=False)

# %%
result_df3[(result_df3['Freq_cancer_IS']-result_df3['Freq_death_notcancer_IS']>=0.1) & (result_df3['integration_variant']=='PhaCinBetalike')]

# %%
genes_t=[]
for gene in result_df3[(result_df3['Freq_cancer_IS']-result_df3['Freq_death_notcancer_IS']>=0.1) & (result_df3['integration_variant']=='KoRV')]['Genes']:
    if gene not in genes_t:
        genes_t.append(gene)
        print(gene) 

# %%
result_df3[(result_df3['Freq_cancer_IS']-result_df3['Freq_death_notcancer_IS']<=-0.1) & (result_df3['integration_variant']=='KoRV')]

# %%
result_df3[(result_df3['Freq_cancer_IS']-result_df3['Freq_death_notcancer_IS']<=-0.1) & (result_df3['integration_variant']=='PhaCinBeta')]

# %%
result_df3[(result_df3['Freq_cancer_IS']-result_df3['Freq_death_notcancer_IS']<=-0.1) & (result_df3['integration_variant']=='PhaCinBetalike')]

# %%
result_df3[(result_df3['Freq_cancer_IS']-result_df3['Freq_death_notcancer_IS']<=-0.1)][["tseq",	"integration_variant","start",
"end","POS","reads_evidence","Missing",	"AAF-alive","AAF-dead","Genes","Freq_cancer_IS","Freq_death_notcancer_IS"]].to_excel(path+"sup_table2.xlsx",index=False)

# %%
genes_t=[]
for gene in result_df3[(result_df3['Freq_cancer_IS']-result_df3['Freq_death_notcancer_IS']<=-0.1) & (result_df3['integration_variant']=='KoRV')]['Genes']:
    if gene not in genes_t:
        genes_t.append(gene)
        print(gene)



