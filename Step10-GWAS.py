# %%
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.utils import resample
import glob
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')
sns.set_theme(style="ticks", palette="colorblind")

# %%
import json
with open('path to GCA_030178435.1/contigs_to_chromosomes.json', 'r') as file:
    chromosomes = json.load(file)
with open('path to GCA_030178435.1/contigs_per_chromosome.json', 'r') as file:
    chro_contigs = json.load(file)
with open('path to GCA_030178435.1/contigs_per_chromosome_postion.json', 'r') as file:
    chro_pos = json.load(file)

# %%
-np.log10(0.05)

# %%
-np.log10(0.05/610867)

# %%
#get annotations
independent_tests = 953591

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
vcf_df[vcf_df["Gene"]=="BCL2L1"]

# %%
#check top results from Paula
threshold_suggestive = -np.log10(0.05)
path = 'path to data'
result_pk = pd.read_csv(path+'R2_plink2.Leukemia.glm.logistic.hybrid',sep='\t')
result_pk['log10_pvalue']=-np.log10(result_pk['P'])
result_pk[(result_pk['log10_pvalue']>=threshold_suggestive) & (result_pk['ID'].str.contains('KoRV'))]

# %%
result_pk[(result_pk['log10_pvalue']>=threshold)]

# %%
result_pk[(result_pk['P']<0.055) & (~result_pk['ID'].str.contains('SNP'))]

# %%
sample_info = pd.read_excel(path + 'Manifest.xlsx')
samples = sample_info["Sample ID"].unique()
sample_info=sample_info.drop_duplicates("Sample ID",keep="first")
pheno = pd.read_excel(path"Koala health outcome data.xlsx",skiprows=1)
len(samples)

# %%
#calculate risk score per individual for offspring
data_leu = result_pk[(result_pk['P']<=0.055) & (~result_pk['ID'].str.contains('SNP'))]
data_leu["ID"]=data_leu["ID"].apply(lambda x: x.split("_")[0]+"_"+x.split("_")[1]+x.split("_")[2])
data_leu['log_OR'] = np.log(data_leu['OR'])

prs_df=pd.merge(data_leu.reset_index(drop=True),vcf_df.reset_index(drop=True),on='ID',how='left')
prs_df=prs_df[~prs_df["Gene"].isna()].reset_index(drop=True)
#count alleles
for i,a1 in enumerate(prs_df['A1']):
    if a1 == prs_df.loc[i,'ALT_x']:
        prs_df.loc[i,:]=prs_df.loc[i,:].replace({'./.':0,'1/1':2,'0/0':0,'1/0':1,'0/1':1})
    else:
        prs_df.loc[i,:]=prs_df.loc[i,:].replace({'./.':0,'1/1':0,'0/0':2,'1/0':1,'0/1':1})

# Calculate PRS for each individual
grs=[]
for sample in sample_info["Sample ID"]:
    grs.append((prs_df['log_OR'] * prs_df[str(sample)]).sum())
grs_df=pd.DataFrame()
grs_df['Accession']=sample_info['Sample ID']
grs_df=pd.merge(grs_df,pheno,how="left",on="Accession")
grs_df['GRS_Leukemia']=grs
grs_df["Leukemia"]=grs_df['COD3'].apply(lambda x: "True" if x==41 else "False")
grs_df

# %%
prs_df

# %%
grs_df.describe()

# %%
# Draw a categorical scatterplot to show each observation
fig, ax = plt.subplots(figsize=(5, 2))
ax = sns.swarmplot(data=grs_df[grs_df["Dead"]==1], x="GRS_Leukemia", y="Sex_c", hue="Leukemia", palette=['darkred','gray'])
ax.set_xlabel('Genetic Risk Score for Leukemia')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, title="Leukemia")
ax.axvline(7, color=sns.color_palette("hls", 8)[0], linestyle='--')
fig.tight_layout()
plt.savefig('grs_leukemia.jpg', dpi=1200)

# %%
# Create a binary column for predictions based on the GRS cutoff of -8
for cutoff in [4,5,6,7,8,9]:
    grs_df['Predicted'] = (grs_df['GRS'] <= cutoff).astype(int)

    # Define number of bootstrap iterations and result storage
    n_iterations = 100
    accuracy_scores = []

    for _ in range(n_iterations):
        # Resample with replacement
        bootstrapped_df = resample(grs_df, replace=True)
        
        # Define features (you can include more features if needed)
        X = bootstrapped_df[['GRS']]  # Independent variable
        y_true = bootstrapped_df['Leukemia'].apply(lambda x: 1 if x=="True" else 0)  # Actual true values based on Offspring column

        # Split the data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(X, y_true, test_size=0.25, random_state=42)

        # Initialize and fit the logistic regression model
        model = LogisticRegression()
        model.fit(X_train, y_train)

        # Make predictions on the test set
        predictions = model.predict(X_test)

        # Calculate accuracy and add to the list
        accuracy = accuracy_score(y_test, predictions)
        accuracy_scores.append(accuracy)

    # Calculate the average accuracy from bootstrapping
    average_accuracy = np.mean(accuracy_scores)
    print(f'Average Prediction Accuracy cutoff {cutoff} bootstraps: {average_accuracy:.4f}')


# %%
grs_df[(grs_df["Dead"]==0) & (grs_df["GRS"]>=7)]

# %%
temp=result_pk[(result_pk['P']<=0.1) & (~result_pk['ID'].str.contains('SNP'))]
#vcf_df=vcf_df.drop_duplicates('ID',keep='first')
temp["ID"]=temp["ID"].apply(lambda x: x.split("_")[0]+"_"+x.split("_")[1]+x.split("_")[2])
temp=pd.merge(temp,vcf_df[['ID','INFO','Impact','Consequence','Gene']],how='left',on='ID')
temp

# %%
samples_info.groupby('Off_total').count()['Accession']

# %%
result = pd.read_csv(path+'R2_plink2.AgeDead.glm.linear',sep='\t')
result=result[~result['P'].isna()]
result["Chromosome"]=result['#CHROM'].replace(chromosomes)

custom_order =[]
for chro in range(1,8):
    temp=pd.DataFrame()
    temp['contig']=chro_contigs[str(chro)]
    temp['pos']=chro_pos[str(chro)]
    temp=temp.astype({'pos':int}).sort_values('pos').reset_index(drop=True)
    custom_order+=list(temp['contig']) # add contigs in ascending order based on mathcing blast positions per chromosome
for chro,contig in chromosomes.items():
    if chro=='X': #add x contigs
        custom_order.append(contig)
for contig in result['#CHROM'].unique():
    if contig not in custom_order: #adding remaining contigs
        custom_order.append(contig)
custom_order_dict = {k: v for v, k in enumerate(custom_order)}
result['CustomOrder'] = result["#CHROM"].map(custom_order_dict)
result=result.astype({"Chromosome":str}).sort_values(by=['CustomOrder']).reset_index(drop=True)
result['ind']=result.index


# Calculate genome-wide significance threshold
threshold = -np.log10(0.05/independent_tests)
result['log10_pvalue']=-np.log10(result['P'])

threshold_sugg = 6

threshold_suggestive = -np.log10(0.05)

fig, ax = plt.subplots(figsize=(17, 3))

x_labels = []
x_labels_pos = []   
colors = {"1":"lightgray","2":"darkgray","3":"lightgray","4":"darkgray","5":"lightgray","6":"darkgray","7":"lightgray","X":"darkgray","Others":"lightgray"}

for chromosome, data in result.groupby('Chromosome'):
    #ax.scatter(data['ind'], data['log10_pvalue'], color="white", alpha=0, s=0)
    data = data[data["ID"].str.contains("SNP")]
    if "JA" in chromosome:
        chromosome="Others"
        label=''
    elif chromosome=='1':
        label='SNP_indel'
    else:
        label=''
    ax.scatter(data['ind'], data['log10_pvalue'], color=colors[chromosome], alpha=0.6, s=3, label=label)
    x_labels.append(chromosome)
    x_labels_pos.append((data['ind'].iloc[-1] - (data['ind'].iloc[-1] - data['ind'].iloc[0])/2))  


data = result[result["ID"].str.contains("KoRV")]
ax.scatter(data['ind'], data['log10_pvalue'], color='indianred',label="KoRV", alpha=1, s=5, marker="D")
    
data = result[result["ID"].str.contains("PhaCinBeta")]
ax.scatter(data['ind'], data['log10_pvalue'], color=sns.color_palette("husl", 8)[1],label=r'phaCin-$\beta$', alpha=1, s=5, marker="D")

data = result[result["log10_pvalue"]>=threshold]
ax.scatter(data['ind'], data['log10_pvalue'], color='red', label='SNP significant', s=10)


data = result[result["log10_pvalue"]>=threshold_sugg]
ax.scatter(data['ind'], data['log10_pvalue'], color='blue', label='SNP suggestive', s=10)

data = result[result["ID"].str.contains("PhaCinBeta") & (result["log10_pvalue"]>=threshold_suggestive)]
ax.scatter(data['ind'], data['log10_pvalue'], color=sns.color_palette("husl", 8)[1],label=r'phaCin-$\beta$'+' suggestive', 
alpha=1, s=20, marker="D",edgecolor='dimgray',linewidths=0.5)

data = result[result["ID"].str.contains("KoRV") & (result["log10_pvalue"]>=threshold_suggestive)]
ax.scatter(data['ind'], data['log10_pvalue'], color='indianred',label="KoRV suggestive", alpha=1, s=20, marker="D",edgecolor='dimgray',linewidths=0.5)

ax.axhline(threshold, color='red', linestyle='--')
ax.axhline(threshold_sugg, color='blue', linestyle='--')
ax.axhline(threshold_suggestive, color='dimgray', linestyle='--')
ax.set(title='Linear model for age at death', xlabel='Virtually Assigned Chromosome', ylabel='$-log_{10}(pvalue)$')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.subplots_adjust(right=0.75)
ax.set_xticks(x_labels_pos[0:8])
ax.set_xticklabels(['1','2','3','4','5','6','7','Others'])

for label in ax.get_xticklabels()[-1:]:#adjust the last label
    label.set_horizontalalignment('left')

ax.set_xlim([0, len(result)])
ax.set_ylim([0, 8])
sns.despine()
plt.savefig('GWAS_longevity.jpg', dpi=1200)
plt.show()

# %%
data_long = result[(result["log10_pvalue"]>=threshold_sugg) | (result["ID"].str.contains("KoRV") & (result["log10_pvalue"]>=threshold_suggestive)) | (result["ID"].str.contains("PhaCinBeta") & (result["log10_pvalue"]>=threshold_suggestive))]
data_long["ID"]=data_long["ID"].apply(lambda x: x.split("_")[0]+"_"+x.split("_")[1]+x.split("_")[2])
data_long =pd.merge(data_long ,vcf_df[['ID','INFO','Impact','Consequence','Gene']],how='left',on='ID')
data_long

# %%
data_long.columns

# %%
len(data_long[data_long["ID"].str.contains("SNP")])

# %%
len(data_long[data_long["ID"].str.contains("KoRV")])

# %%
len(data_long[data_long["ID"].str.contains("PhaCinBeta")])

# %%
snps_replace = {"US_"+str(sample):str(sample) for sample in samples}
snps.rename(snps_replace,axis=1).rename({"ID_x":"ID","ALT_x":"ALT"},axis=1)

# %%
data_long2 =pd.merge(data_long,vcf_df[['ID']+[str(sample) for sample in samples]],how='left',on='ID')
snps_replace = {"US_"+str(sample):str(sample) for sample in samples}
data_long2 = pd.concat([data_long2[~data_long2["ID"].str.contains("SNP")],snps.rename(snps_replace,axis=1).rename({"ID_x":"ID","ALT_x":"ALT"},axis=1)])
data_long2[~data_long2["Gene"].isna()]

# %%
prs_df = data_long2[~data_long2["Gene"].isna()].reset_index(drop=True)
len(prs_df[prs_df["ID"].str.contains("KoRV")])

# %%
len(prs_df[prs_df["ID"].str.contains("SNP")])

# %%
len(prs_df[prs_df["ID"].str.contains("PhaCinBeta")])

# %%
prs_df = data_long2[~data_long2["Gene"].isna()].reset_index(drop=True)

# Count alleles
for i, a1 in enumerate(prs_df['A1']):
    if a1 == prs_df.loc[i, 'ALT']:
        prs_df.loc[i, :] = prs_df.loc[i, :].replace({'./.': 0, '1/1': 2, '0/0': 0, '1/0': 1, '0/1': 1})
    else:
        prs_df.loc[i, :] = prs_df.loc[i, :].replace({'./.': 0, '1/1': 0, '0/0': 2, '1/0': 1, '0/1': 1})

# Calculate FRS (Family Risk Score) for each individual
# The FRS is calculated as the sum of beta values multiplied by the count of alleles
frs = []
for sample in sample_info['Sample ID']:
    frs.append((prs_df['BETA'] * prs_df[str(sample)]).sum())  # Use beta instead of log_OR

# Create a DataFrame to store the results
frs_df = pd.DataFrame()
frs_df['Accession']=sample_info['Sample ID']
frs_df=pd.merge(frs_df,pheno,how="left",on="Accession")
frs_df['FRS'] = frs
frs_df

# %%
frs_df[frs_df["Dead"]==1].describe()

# %%

from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_curve, auc

# Assuming frs_df includes a column 'Trait' for the continuous outcome
X =frs_df[frs_df["Dead"]==1][['FRS']]
y = frs_df[frs_df["Dead"]==1]['Age1']  # Your continuous trait of interest

# Split the data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42)

# Fit a simple linear regression model
model = LinearRegression()
model.fit(X_train, y_train)

# Predicted values
y_pred = model.predict(X_test)

# Determine an arbitrary threshold using y_pred for binary classification, e.g., median split
threshold = np.median(y_pred)
y_pred_binary = (y_pred >= threshold).astype(int)
y_test_median_split = (y_test >= np.median(y_test)).astype(int)

# Compute ROC curve and AUC
fpr, tpr, _ = roc_curve(y_test_median_split, y_pred_binary)
roc_auc = auc(fpr, tpr)

# Plot ROC Curve
plt.figure()
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.show()


# %%
frs_df["Sex"]=frs_df["Sex_c"]
g = sns.lmplot(
    data=frs_df[frs_df["Dead"] == 1], 
    x="Age1", 
    y="FRS", 
    hue="Sex",
    palette="muted", 
    ci=None, 
    height=4, 
    scatter_kws={"s": 50, "alpha": 1}
)

# Use the FacetGrid methods to set labels
g.set_axis_labels('Age at Death', 'Index')
g.set_titles('Breeding Index for Longevity')

# If you need to adjust the layout or appearance further
plt.subplots_adjust(top=0.9)
plt.show()
plt.savefig('index_longevity.jpg', dpi=1200)


# %%
vcf_file_path = path+'DBimportgeno_consensus.annotated.vcf' 
vcf_snps = read_vcf(vcf_file_path)
vcf_snps['Gene']=vcf_snps["INFO"].apply(lambda x: x.split("|")[3])
vcf_snps['Impact']=vcf_snps["INFO"].apply(lambda x: x.split("|")[2])
vcf_snps['Consequence']=vcf_snps["INFO"].apply(lambda x: x.split("|")[1])
vcf_snps.head()

# %%
snps=data_long[data_long["ID"].str.contains("SNP")].sort_values("POS")
snps=snps.astype({"#CHROM":str,"POS":str}).drop(["Gene","Impact","Consequence"],axis=1)
snps=pd.merge(snps,vcf_snps.astype({"#CHROM":str,"POS":str}),how="left",on=["#CHROM","POS"])
snps

# %%
result = pd.read_csv(path+'R2_plink2.OffspringTotalBinDead.glm.logistic.hybrid',sep='\t')
result=result[~result['P'].isna()]
result["Chromosome"]=result['#CHROM'].replace(chromosomes)

custom_order =[]
for chro in range(1,8):
    temp=pd.DataFrame()
    temp['contig']=chro_contigs[str(chro)]
    temp['pos']=chro_pos[str(chro)]
    temp=temp.astype({'pos':int}).sort_values('pos').reset_index(drop=True)
    custom_order+=list(temp['contig']) # add contigs in ascending order based on mathcing blast positions per chromosome
for chro,contig in chromosomes.items():
    if chro=='X': #add x contigs
        custom_order.append(contig)
for contig in result['#CHROM'].unique():
    if contig not in custom_order: #adding remaining contigs
        custom_order.append(contig)
custom_order_dict = {k: v for v, k in enumerate(custom_order)}
result['CustomOrder'] = result["#CHROM"].map(custom_order_dict)
result=result.astype({"Chromosome":str}).sort_values(by=['CustomOrder']).reset_index(drop=True)
result['ind']=result.index


# Calculate genome-wide significance threshold
threshold = -np.log10(0.05/independent_tests)
result['log10_pvalue']=-np.log10(result['P'])

threshold_sugg = 6

threshold_suggestive = -np.log10(0.05)

fig, ax = plt.subplots(figsize=(17, 3))

x_labels = []
x_labels_pos = []   
colors = {"1":"lightgray","2":"darkgray","3":"lightgray","4":"darkgray","5":"lightgray","6":"darkgray","7":"lightgray","X":"darkgray","Others":"lightgray"}

for chromosome, data in result.groupby('Chromosome'):
    #ax.scatter(data['ind'], data['log10_pvalue'], color="white", alpha=0, s=0)
    data = data[data["ID"].str.contains("SNP")]
    if "JA" in chromosome:
        chromosome="Others"
        label=''
    elif chromosome=='1':
        label='SNP_indel'
    else:
        label=''
    ax.scatter(data['ind'], data['log10_pvalue'], color=colors[chromosome], alpha=0.6, s=3, label=label)
    x_labels.append(chromosome)
    x_labels_pos.append((data['ind'].iloc[-1] - (data['ind'].iloc[-1] - data['ind'].iloc[0])/2))  


data = result[result["ID"].str.contains("KoRV")]
ax.scatter(data['ind'], data['log10_pvalue'], color='indianred',label="KoRV", alpha=1, s=5, marker="D")
    
data = result[result["ID"].str.contains("PhaCinBeta")]
ax.scatter(data['ind'], data['log10_pvalue'], color=sns.color_palette("husl", 8)[1],label=r'phaCin-$\beta$', alpha=1, s=5, marker="D")

data = result[result["log10_pvalue"]>=threshold]
ax.scatter(data['ind'], data['log10_pvalue'], color='red', label='SNP significant', s=10)


data = result[result["log10_pvalue"]>=threshold_sugg]
ax.scatter(data['ind'], data['log10_pvalue'], color='blue', label='SNP suggestive', s=10)

data = result[result["ID"].str.contains("PhaCinBeta") & (result["log10_pvalue"]>=threshold_suggestive)]
ax.scatter(data['ind'], data['log10_pvalue'], color=sns.color_palette("husl", 8)[1],label=r'phaCin-$\beta$'+' suggestive', 
alpha=1, s=20, marker="D",edgecolor='dimgray',linewidths=0.5)

data = result[result["ID"].str.contains("KoRV") & (result["log10_pvalue"]>=threshold_suggestive)]
ax.scatter(data['ind'], data['log10_pvalue'], color='indianred',label="KoRV suggestive", alpha=1, s=20, marker="D",edgecolor='dimgray',linewidths=0.5)

ax.axhline(threshold, color='red', linestyle='--')
ax.axhline(threshold_sugg, color='blue', linestyle='--')
ax.axhline(threshold_suggestive, color='dimgray', linestyle='--')
ax.set(title='Logistic model for having an offspring', xlabel='Virtually Assigned Chromosome', ylabel='$-log_{10}(pvalue)$')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.subplots_adjust(right=0.75)
ax.set_xticks(x_labels_pos[0:8])
ax.set_xticklabels(['1','2','3','4','5','6','7','Others'])

for label in ax.get_xticklabels()[-1:]:#adjust the last label
    label.set_horizontalalignment('left')

ax.set_xlim([0, len(result)])
ax.set_ylim([0, 8])
sns.despine()
plt.savefig('GWAS_fertility_PK.jpg', dpi=1200)
plt.show()

# %%
result = pd.read_csv(path+'R2_plink2.OffspringTotalBinDead.glm.logistic.hybrid',sep='\t')
result=result[~result['P'].isna()]
result["Chromosome"]=result['#CHROM'].replace(chromosomes)

custom_order =[]
for chro in range(1,8):
    temp=pd.DataFrame()
    temp['contig']=chro_contigs[str(chro)]
    temp['pos']=chro_pos[str(chro)]
    temp=temp.astype({'pos':int}).sort_values('pos').reset_index(drop=True)
    custom_order+=list(temp['contig']) # add contigs in ascending order based on mathcing blast positions per chromosome
for chro,contig in chromosomes.items():
    if chro=='X': #add x contigs
        custom_order.append(contig)
for contig in result['#CHROM'].unique():
    if contig not in custom_order: #adding remaining contigs
        custom_order.append(contig)
custom_order_dict = {k: v for v, k in enumerate(custom_order)}
result['CustomOrder'] = result["#CHROM"].map(custom_order_dict)
result=result.astype({"Chromosome":str}).sort_values(by=['CustomOrder']).reset_index(drop=True)
result['ind']=result.index


# Calculate genome-wide significance threshold
threshold = -np.log10(0.05/independent_tests)
result['log10_pvalue']=-np.log10(result['P'])

threshold_sugg = 6

threshold_suggestive = -np.log10(0.05)

# %%
data_off = result[(result["log10_pvalue"]>=threshold_sugg) | (result["ID"].str.contains("KoRV") & (result["log10_pvalue"]>=threshold_suggestive)) | (result["ID"].str.contains("PhaCinBeta") & (result["log10_pvalue"]>=threshold_suggestive))]
data_off["ID"]=data_off["ID"].apply(lambda x: x.split("_")[0]+"_"+x.split("_")[1]+x.split("_")[2])
data_off =pd.merge(data_off ,vcf_df[['ID','INFO','Impact','Consequence','Gene']],how='left',on='ID')
data_off

# %%
#calculate risk score per individual for offspring
data_off['log_OR'] = np.log(data_off['OR'])

prs_df=pd.merge(data_off.reset_index(drop=True),vcf_df.reset_index(drop=True),on='ID',how='left')

#count alleles
for i,a1 in enumerate(prs_df['A1']):
    if a1 == prs_df.loc[i,'ALT_x']:
        prs_df.loc[i,:]=prs_df.loc[i,:].replace({'./.':0,'1/1':2,'0/0':0,'1/0':1,'0/1':1})
    else:
        prs_df.loc[i,:]=prs_df.loc[i,:].replace({'./.':0,'1/1':0,'0/0':2,'1/0':1,'0/1':1})

# Calculate PRS for each individual
grs=[]
for sample in sample_info["Sample ID"]:
    grs.append((prs_df['log_OR'] * prs_df[str(sample)]).sum())
grs_df2=pd.DataFrame()
grs_df2['Accession']=sample_info['Sample ID']
grs_df2=pd.merge(grs_df2,pheno,how="left",on="Accession")
grs_df2['GRS']=grs
grs_df2["Offspring"]=grs_df2['Off_total'].apply(lambda x: "True" if int(x)>0 else "False")
grs_df2

# %%
grs_df['GRS_Offspring']=grs
grs_df["Offspring"]=grs_df2["Offspring"]
#grs_df["Breeding_index_longevity"]=frs
grs_df

# %%
grs_df.to_excel("GRS.xlsx",index=False)

# %%
# Draw a categorical scatterplot to show each observation
fig, ax = plt.subplots(figsize=(5, 2))
ax = sns.swarmplot(data=grs_df[grs_df["Dead"]==1], x="GRS_Offspring", y="Sex_c", hue="Offspring", palette=['gray','green'])
ax.set_xlabel('Genetic Risk Score for Reproduction Success')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, title="Reproduction \nSuccess")
ax.axvline(-7, color=sns.color_palette("hls", 8)[0], linestyle='--')
fig.tight_layout()
plt.savefig('grs_offspring.jpg', dpi=1200)

# %%
grs_df[grs_df["Dead"]==0].describe()

# %%


# Assuming grs_df is your DataFrame with GRS already calculated and the Offspring column exists

# Create a binary column for predictions based on the GRS cutoff of -8
for cutoff in [-10, -9, -8, -7, -6, -5]:
    grs_df['Predicted'] = (grs_df['GRS'] <= cutoff).astype(int)

    # Define number of bootstrap iterations and result storage
    n_iterations = 100
    accuracy_scores = []

    for _ in range(n_iterations):
        # Resample with replacement
        bootstrapped_df = resample(grs_df, replace=True)
        
        # Define features (you can include more features if needed)
        X = bootstrapped_df[['GRS']]  # Independent variable
        y_true = bootstrapped_df['Offspring'].apply(lambda x: 1 if x=="True" else 0)  # Actual true values based on Offspring column

        # Split the data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(X, y_true, test_size=0.25, random_state=42)

        # Initialize and fit the logistic regression model
        model = LogisticRegression()
        model.fit(X_train, y_train)

        # Make predictions on the test set
        predictions = model.predict(X_test)

        # Calculate accuracy and add to the list
        accuracy = accuracy_score(y_test, predictions)
        accuracy_scores.append(accuracy)

    # Calculate the average accuracy from bootstrapping
    average_accuracy = np.mean(accuracy_scores)
    print(f'Average Prediction Accuracy cutoff {cutoff} bootstraps: {average_accuracy:.4f}')


# %%
grs_df[(grs_df["Dead"]==0) & (grs_df["GRS"]<=-7)][["Accession","Name","GRS","Age1","Sex_c","Offspring","Off_total"]]

# %%
samples_info['Offspring']=samples_info['Off_total'].apply(lambda x: 1 if x==0 else 2) #24 no offsrping, 34 with some offspring
samples_info[['Accession','Sex','Cancer','Age1','Offspring']].rename(columns={'Sex':'sex','Age1':'Age',
    'Accession': '#IID'}).to_csv('pheno_fertility.txt',sep='\t',index=False)
samples_info.groupby('Offspring').count()['Accession']

# %%
# Draw a categorical scatterplot to show each observation
fig, ax = plt.subplots(figsize=(5, 2))
ax = sns.swarmplot(data=grs_df[grs_df["Dead"]==0], x="GRS_Leukemia", y="Sex_c", color='gray')
ax.set_xlabel('Genetic Risk Score for Leukemia')
ax.axvline(7, color=sns.color_palette("hls", 8)[0], linestyle='--')
fig.tight_layout()
plt.savefig('grs_leukemia_living.jpg', dpi=1200)

# %%
fig, ax = plt.subplots(figsize=(5, 2))
ax = sns.swarmplot(data=grs_df[grs_df["Dead"]==0], x="GRS_Offspring", y="Sex_c", hue="Offspring", palette=[sns.color_palette("hls", 8)[0],'gray'])
ax.set_xlabel('Genetic Risk Score for Fertility')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, title="Offspring")
ax.axvline(-7, color=sns.color_palette("hls", 8)[0], linestyle='--')
fig.tight_layout()
plt.savefig('grs_offspring_living.jpg', dpi=1200)


