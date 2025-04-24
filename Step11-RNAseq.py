# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %%
data = pd.read_csv("path to DE_analysis.csv").rename({"Unnamed: 0":"Gene"},axis=1)
data

# %%
data["pvalue"].isna().sum()

# %%
data["padj"].isna().sum()

# %%
data[data['baseMean']>=10]["padj"].isna().sum()

# %%
annot = pd.read_csv("path to GCA_030178435.1_masked_KoRV_PhER_PhaCinBeta/braker/braker.gtf",sep="\t",
names=['Contig','Source','Type','start','end','col','strand','col2','INFO'])
contigs={}
contigs_f = open("GCA_030178435.1_ASM3017843v1_genomic.fna.ann","r")
for line in contigs_f.readlines():
    if "whole genome shotgun sequence" in line:
        contig=line.split(" ")[6][:-1]
        name=line.split(" ")[1]
        contigs[contig]=name
annot['Contig']=annot['Contig'].replace(contigs)
contigs_f.close()
orthologs = pd.read_excel("path to GCA_030178435.1/Supp_Table2.xlsx",comment="#")
orthologs = orthologs[~orthologs['query'].isna()]
orthologs['query']=orthologs['query'].apply(lambda x: x.split(".")[0])
genes = annot[annot['Type']=='gene'].reset_index(drop=True)
genes = pd.merge(orthologs,genes,how='left',left_on='query',right_on='INFO')
genes

# %%
data=data[(~data["padj"].isna()) & (data['baseMean']>=10)]
genes.index=genes["query"]
data['-log10_pvalue'] = -np.log10(data['padj'])
data['gene']=data["Gene"].apply(lambda x: genes.loc[x.split(".")[0],"Preferred_name"] if x.split(".")[0] in genes.index else x)
# Identify DEGs
# Adjust the threshold according to your specific criteria
log2fc_threshold = 1.0
pvalue_threshold = 0.05

data['is_DEG'] = (((data['log2FoldChange'] >= log2fc_threshold) & (data['padj'] <= pvalue_threshold)) |
((data['log2FoldChange'] <= -log2fc_threshold) & (data['padj'] <= pvalue_threshold)) )

# Plot
plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=data,
    x='log2FoldChange',
    y='-log10_pvalue',
    hue='is_DEG',
    palette={True: "indianred", False: "grey"},
    legend=False
)

# Label the DEGs
for _, row in data.iterrows():
    if row['is_DEG'] and (row['-log10_pvalue']>=25 or row["log2FoldChange"]>13) and type(row['gene'])==str:
        plt.text(row['log2FoldChange'], row['-log10_pvalue'], row['gene'], horizontalalignment='left',fontsize=8,rotation=30)

# Add labels and title
plt.axhline(y=-np.log10(pvalue_threshold), color='black', linestyle='--', linewidth=0.8)
plt.axvline(x=log2fc_threshold, color='black', linestyle='--', linewidth=0.8)
plt.axvline(x=-log2fc_threshold, color='black', linestyle='--', linewidth=0.8)
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 P-value')
plt.title('Volcano Plot of Differential Expression')

plt.show()


# %%
eu_is=pd.read_excel(path+"IS_detected_annoattedEU.xlsx")
eu_is=eu_is[(eu_is["395"]==1) |(eu_is["463"]==1) | (eu_is["407"]==1) | (eu_is["436"]==1) | (eu_is["461"]==1)| (eu_is["399"]==1) ]
eu_is

# %%
korv_genes=[]
for gene in eu_is[~eu_is["Genes"].isna()]["Genes"]:
    for item in gene.split(","):
        if item.strip() not in korv_genes:
            korv_genes.append(item.strip())
len(korv_genes)

# %%
eu_is['count']=eu_is["395"]+eu_is["463"]+eu_is["407"]+eu_is["436"]+eu_is["461"]+eu_is["399"]
korv_genes_high=[]
for gene in eu_is[(~eu_is["Genes"].isna()) & (eu_is['count']>=2)]["Genes"]:
    for item in gene.split(","):
        if item.strip() not in korv_genes_high:
            korv_genes_high.append(item.strip())
len(korv_genes_high)

# %%
korv_genes_high

# %%
data['is_DEG'].sum()

# %%
len(data)

# %%
plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=data,
    x='log2FoldChange',
    y='-log10_pvalue',
    hue='is_DEG',
    palette={True: "indianred", False: "grey"},
    legend=False
)

# Label the DEGs
for _, row in data.iterrows():
    if type(row['gene'])==str:
        if row['is_DEG'] and (row['-log10_pvalue']>=30 or row["log2FoldChange"]<-12):
            plt.text(row['log2FoldChange'], row['-log10_pvalue'], row['gene'], rotation=30,fontsize=8, color='black')
        if row['is_DEG'] and (row['gene'] in korv_genes_high):
            print(row['gene'])
            plt.text(row['log2FoldChange'], row['-log10_pvalue'], row['gene'], horizontalalignment='left', fontsize=8, color='black', fontweight='bold',
         bbox=dict(facecolor='yellow', alpha=0.5))

# Add labels and title
plt.axhline(y=-np.log10(pvalue_threshold), color='black', linestyle='--', linewidth=0.8)
plt.axvline(x=log2fc_threshold, color='black', linestyle='--', linewidth=0.8)
plt.axvline(x=-log2fc_threshold, color='black', linestyle='--', linewidth=0.8)
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 P-value')
plt.title('Volcano Plot of Differential Expression')

plt.show()

# %%
plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=data,
    x='log2FoldChange',
    y='-log10_pvalue',
    hue='is_DEG',
    palette={True: "indianred", False: "grey"},
    legend=False
)

# Label the DEGs
for _, row in data.iterrows():
    if type(row['gene'])==str:
        if row['is_DEG'] and (row['-log10_pvalue']>=30 or row["log2FoldChange"]>12):
            plt.text(row['log2FoldChange'], row['-log10_pvalue'], row['gene'], rotation=30,fontsize=8, color='black')
        if row['is_DEG'] and (row['gene'] in korv_genes):
            print(row['gene'])
            plt.text(row['log2FoldChange'], row['-log10_pvalue'], row['gene'], horizontalalignment='right',fontsize=8, color='black', fontweight='bold',
         bbox=dict(facecolor='yellow', alpha=0.5))

# Add labels and title
plt.axhline(y=-np.log10(pvalue_threshold), color='black', linestyle='--', linewidth=0.8)
plt.axvline(x=log2fc_threshold, color='black', linestyle='--', linewidth=0.8)
plt.axvline(x=-log2fc_threshold, color='black', linestyle='--', linewidth=0.8)
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 P-value')
plt.title('Volcano Plot of Differential Expression')

plt.show()

# %%
data["gene"]=data["gene"].fillna("")
data[(data["gene"].str.contains("SLC29A1", na=False))]

# %%
data[(data["gene"].str.contains("GAPDH", na=False))]

# %%
genes.loc["g15951",:]

# %%
genes.loc["g15952",:]

# %%
counts = pd.read_csv(path+"Normalized_counts.csv").rename({"Unnamed: 0":"Gene"},axis=1)
counts['gene']=counts["Gene"].apply(lambda x: genes.loc[x.split(".")[0],"Preferred_name"] if x.split(".")[0] in genes.index else x)
counts

# %%
counts[(counts["gene"].str.contains("SLC29A1", na=False))]

# %%
eu_is[eu_is['Genes']=='SLC29A1']

# %%
from scipy.stats import mannwhitneyu
temp=pd.DataFrame(counts.loc[22650,counts.columns[1:-1]].T)
status=[]
for i in temp.index:
    if 'ER' not in i and eu_is.loc[1075,i.split('_')[0]]==1:
        status.append('Integration')
    else:
        status.append('No Integration')

temp['KoRV Status']=status
integration_counts = list(temp[temp['KoRV Status'] == 'Integration'][22650])
no_integration_counts = list(temp[temp['KoRV Status'] == 'No Integration'][22650])

u_stat, p_value = mannwhitneyu(integration_counts, no_integration_counts, alternative='two-sided')

sns.boxplot(x="KoRV Status", y=22650,data=temp,palette=['mediumaquamarine','indianred'], fliersize=2,linewidth=3)
sns.stripplot(temp, x="KoRV Status", y=22650, size=7, color=".3",jitter=False,linewidth=.5,edgecolor='black')
# Determine the positions for annotation
y, h, col = temp[22650].max() + 50, 5, 'k'

# Annotate p-value on the plot
plt.text(0.5, y + h, f'p = {p_value:.4f}', ha='center', va='bottom', fontsize=12, color='black')

# Draw a line below the p-value annotation
plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1, color='.3')
plt.title("Expression of SLC29A1 Containing an Intronic KoRV Integration")
plt.xlabel("")

plt.ylabel("Normalized Count")
sns.despine(offset=10, trim=True)

# %%
counts[(counts["gene"].str.contains("PCCA", na=False))]

# %%
eu_is[eu_is['Genes']=='PCCA']

# %%
from scipy.stats import mannwhitneyu
temp=pd.DataFrame(counts.loc[9616,counts.columns[1:-1]].T)
status=[]
for i in temp.index:
    if 'ER' not in i and eu_is.loc[1512,i.split('_')[0]]==1:
        status.append('Integration')
    else:
        status.append('No Integration')

temp['KoRV Status']=status
integration_counts = list(temp[temp['KoRV Status'] == 'Integration'][9616])
no_integration_counts = list(temp[temp['KoRV Status'] == 'No Integration'][9616])

u_stat, p_value = mannwhitneyu(integration_counts, no_integration_counts, alternative='two-sided')

sns.boxplot(x="KoRV Status", y=9616,data=temp,palette=['mediumaquamarine','indianred'], fliersize=2,linewidth=3)
sns.stripplot(temp, x="KoRV Status", y=9616, size=7, color=".3",jitter=False,linewidth=.5,edgecolor='black')
# Determine the positions for annotation
y, h, col = temp[9616].max() + 50, 5, 'k'

# Annotate p-value on the plot
plt.text(0.5, y + h, f'p = {p_value:.4f}', ha='center', va='bottom', fontsize=12, color='black')

# Draw a line below the p-value annotation
plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1, color='.3')
plt.title("Expression of PCCA Containing a KoRV Integration")
plt.xlabel("")

plt.ylabel("Normalized Count")
sns.despine(offset=10, trim=True)

# %%
counts[(counts["gene"].str.contains("TRAPPC12", na=False))]

# %%
eu_is[eu_is['Genes']=='TRAPPC12']

# %%
col=20542
temp=pd.DataFrame(counts.loc[col,counts.columns[1:-1]].T)
status=[]
for i in temp.index:
    if 'ER' not in i and eu_is.loc[585,i.split('_')[0]]==1:
        status.append('Integration')
    else:
        status.append('No Integration')

temp['KoRV Status']=status
integration_counts = list(temp[temp['KoRV Status'] == 'Integration'][col])
no_integration_counts = list(temp[temp['KoRV Status'] == 'No Integration'][col])

u_stat, p_value = mannwhitneyu(integration_counts, no_integration_counts, alternative='two-sided')

sns.boxplot(x="KoRV Status", y=col,data=temp,palette=['mediumaquamarine','indianred'], fliersize=2,linewidth=3)
sns.stripplot(temp, x="KoRV Status", y=col, size=7, color=".3",jitter=False,linewidth=.5,edgecolor='black')
# Determine the positions for annotation
y, h, col = temp[col].max() + 50, 5, 'k'

# Annotate p-value on the plot
plt.text(0.5, y + h, f'p = {p_value:.4f}', ha='center', va='bottom', fontsize=12, color='black')

# Draw a line below the p-value annotation
plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1, color='.3')
plt.title("Expression of TRAPPC12 Containing a KoRV Integration")
plt.xlabel("")

plt.ylabel("Normalized Count")
sns.despine(offset=10, trim=True)

# %%
counts[(counts["gene"].str.contains("SERINC5", na=False))]

# %%
eu_is[eu_is['Genes']=='SERINC5']

# %%
col=2261
temp=pd.DataFrame(counts.loc[col,counts.columns[1:-1]].T)
status=[]
for i in temp.index:
    if 'ER' not in i and eu_is.loc[837,i.split('_')[0]]==1:
        status.append('Integration')
    else:
        status.append('No Integration')

temp['KoRV Status']=status
integration_counts = list(temp[temp['KoRV Status'] == 'Integration'][col])
no_integration_counts = list(temp[temp['KoRV Status'] == 'No Integration'][col])

u_stat, p_value = mannwhitneyu(integration_counts, no_integration_counts, alternative='two-sided')

sns.boxplot(x="KoRV Status", y=col,data=temp,palette=['mediumaquamarine','indianred'], fliersize=2,linewidth=3)
sns.stripplot(temp, x="KoRV Status", y=col, size=7, color=".3",jitter=False,linewidth=.5,edgecolor='black')
# Determine the positions for annotation
y, h, col = temp[col].max() + 50, 5, 'k'

# Annotate p-value on the plot
plt.text(0.5, y + h, f'p = {p_value:.4f}', ha='center', va='bottom', fontsize=12, color='black')

# Draw a line below the p-value annotation
plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1, color='.3')
plt.title("Expression of SERINC5 Containing a KoRV Integration")
plt.xlabel("")

plt.ylabel("Normalized Count")
sns.despine(offset=10, trim=True)

# %%
counts[(counts["gene"].str.contains("ING2", na=False))]

# %%
eu_is[eu_is['Genes']=='ING2']

# %%
eu_is[eu_is['Genes']=='LONP2']

# %%
from scipy.stats import mannwhitneyu
col=6544
temp=pd.DataFrame(counts.loc[col,counts.columns[1:-1]].T)
status=[]
for i in temp.index:
    if 'ER' not in i and eu_is.loc[1320,i.split('_')[0]]==1:
        status.append('Integration')
    else:
        status.append('No Integration')

temp['KoRV Status']=status
integration_counts = list(temp[temp['KoRV Status'] == 'Integration'][col])
no_integration_counts = list(temp[temp['KoRV Status'] == 'No Integration'][col])

u_stat, p_value = mannwhitneyu(integration_counts, no_integration_counts, alternative='two-sided')

sns.boxplot(x="KoRV Status", y=col,data=temp,palette=['mediumaquamarine','indianred'], fliersize=2,linewidth=3)
sns.stripplot(temp, x="KoRV Status", y=col, size=7, color=".3",jitter=False,linewidth=.5,edgecolor='black')
# Determine the positions for annotation
y, h, col = temp[col].max() + 50, 5, 'k'

# Annotate p-value on the plot
plt.text(0.5, y + h, f'p = {p_value:.4f}', ha='center', va='bottom', fontsize=12, color='black')

# Draw a line below the p-value annotation
plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1, color='.3')
plt.title("Expression of ING2 Containing a KoRV Integration")
plt.xlabel("")

plt.ylabel("Normalized Count")
sns.despine(offset=10, trim=True)

# %%
plt.rcParams.update({
    'axes.titlesize': 16,     # Title font size
    'axes.labelsize': 16,     # X and Y label font size
    'xtick.labelsize': 16,    # X tick labels font size
    'ytick.labelsize': 16,    # Y tick labels font size
    'legend.fontsize': 12,    # Legend font size
    'font.size': 16           # Base font size for other text
})

for gene in ['PLAT','LONP2','STAP1','SLC29A1','PPM1H','LARP7','PCSK2','SH3KBP1','LARP7','RPL15','ING2','MRPL39','PCCA']:
    try:
        col=list(counts[(counts["gene"].str.contains(gene, na=False))].index)[-1]
        temp=pd.DataFrame(counts.loc[col,counts.columns[1:-1]].T)
        status=[]
        for i in temp.index:
            if 'ER' not in i and eu_is.loc[list(eu_is[eu_is['Genes']==gene].index)[0],i.split('_')[0]]==1:
                status.append('Integration')
            else:
                status.append('No Integration')
    except:
        continue

    temp['KoRV Status']=status
    integration_counts = list(temp[temp['KoRV Status'] == 'Integration'][col])
    no_integration_counts = list(temp[temp['KoRV Status'] == 'No Integration'][col])

    u_stat, p_value = mannwhitneyu(integration_counts, no_integration_counts, alternative='two-sided')
    plt.figure(figsize=(5, 3))
    sns.boxplot(x="KoRV Status", y=col,data=temp,palette=['mediumaquamarine','indianred'], fliersize=2,linewidth=3)
    sns.stripplot(temp, x="KoRV Status", y=col, size=7, color=".3",jitter=False,linewidth=.5,edgecolor='black')
    # Determine the positions for annotation
    y, h, col = temp[col].max() + 50, 5, 'k'

    # Annotate p-value on the plot
    plt.text(0.5, y + h, f'p = {p_value:.4f}', ha='center', va='bottom', fontsize=15, color='black')

    # Draw a line below the p-value annotation
    plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1, color='.3')
    plt.title("Expression of "+gene+" Containing an enKoRV \n")
    plt.xlabel("")

    plt.ylabel("Normalized Count")
    sns.despine(offset=10, trim=True)

# %%
counts[(counts["gene"].str.contains("GAPDH", na=False))]

# %%
genes[genes['Preferred_name'].str.contains('SLC29A1')]

# %%
temp_GAPDH=pd.DataFrame(counts.loc[3258,counts.columns[1:-1]].T)
temp_GAPDH['KoRV Status']=["Positive (n=6)"]*6 + ["Negative (n=5)"]*5
sns.boxplot(x="KoRV Status", y=3258,data=temp_GAPDH)
sns.stripplot(temp_GAPDH, x="KoRV Status", y=3258, size=6, color=".3",jitter=False)
plt.title("GAPDH")
plt.ylabel("Normalized Count")
sns.despine(offset=10, trim=True)

# %%
temp["SLC29A1"]=temp[22650]-temp_GAPDH[3258]
sns.boxplot(x="KoRV Status", y="SLC29A1",data=temp,palette="pastel", width=0.6, showfliers=False)
sns.stripplot(temp, x="KoRV Status", y="SLC29A1", size=6, color=".3",jitter=False)
plt.title("SLC29A1")
plt.ylabel("Normalized Count")
sns.despine(offset=10, trim=True)

# %%
temp=pd.DataFrame(counts.loc[22650,counts.columns[1:-1]].T)
temp['KoRV Status']=["Positive (n=6)"]*6 + ["Negative (n=5)"]*5
sns.boxplot(x="KoRV Status", y=22650,data=temp)
sns.stripplot(temp, x="KoRV Status", y=22650, size=6, color=".3",jitter=False)
plt.title("SLC29A1")
plt.ylabel("Normalized Count")
sns.despine(offset=10, trim=True)

# %%
genes[genes['Preferred_name']=="SLC29A1"]

# %%
counts[(counts["gene"].str.contains("PLAT", na=False))]

# %%
temp=pd.DataFrame(counts.loc[29716,counts.columns[1:-1]].T)
temp['KoRV Status']=["Positive (n=6)"]*6 + ["Negative (n=5)"]*5
sns.boxplot(x="KoRV Status", y=29716,data=temp)
sns.stripplot(temp, x="KoRV Status", y=29716, size=6, color=".3",jitter=False)
plt.title("PLAT")
plt.ylabel("Normalized Count")
sns.despine(offset=10, trim=True)

# %%
counts[(counts["gene"].str.contains("PCSK2", na=False))]

# %%
counts[(counts["gene"].str.contains("PPM1H", na=False))]

# %%
counts[(counts["gene"].str.contains("SH3KBP1", na=False))]

# %%
counts[(counts["gene"].str.contains("FANCI", na=False))]

# %%
counts[(counts["gene"].str.contains("FAM71B", na=False))]


