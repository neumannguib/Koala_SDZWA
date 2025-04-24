# %%
import pandas as pd
import allel
import os
import allel.chunked
allel.chunked.storage_registry['default'] = allel.chunked.hdf5mem_zlib1_storage

# %% [markdown]
# Germline mutation rate 

# %%
#get triads and coverages
pheno = pd.read_excel("Koala health outcome data.xlsx",skiprows=1)

# %%
keep = set()
trios = []
pheno.index=pheno['IndividualId']
for sample in pheno['IndividualId']:
    if (pheno.loc[sample,'Father'] in list(pheno['IndividualId'])) and (pheno.loc[sample,'Mother'] in list(pheno['IndividualId'])) :
        keep.update((sample,pheno.loc[sample,'Father'],pheno.loc[sample,'Mother']))
        trios.append((pheno.loc[sample,'Accession'], pheno.loc[pheno.loc[sample,'Father'],'Accession'] , pheno.loc[pheno.loc[sample,'Mother'],'Accession'], pheno.loc[sample,'Sex'] ))

pheno = pheno.loc[list(keep),:]


# Write the PED file
with open(output_file, "w") as file:
    for i, trio in enumerate(trios, start=1):
        family_id = f"FAM{i}"
        offspring_id, sire_id, dam_id, sex = trio
        # Write trio information to the PED file
        if offspring_id==590427:
            file.write(f"{family_id} US_{offspring_id}_liver US_{sire_id} US_{dam_id} {sex} 0\n")
            file.write(f"{family_id} US_{offspring_id}_lung US_{sire_id} US_{dam_id} {sex} 0\n")
        else:
            file.write(f"{family_id} US_{offspring_id} US_{sire_id} US_{dam_id} {sex} 0\n")

        

# %%
os.system('vcftools --gzvcf all_contigs_filteredv2.recode.vcf.gz \
--out unique_variants --recode --recode-INFO-all --max-missing 0.9 --remove-filtered-all --mac 1 --max-mac 2')

# %%
os.system("bgzip unique_variants.recode.vcf")
os.system("tabix unique_variants.recode.vcf.gz")
#run gatk to count for mendelian incosistencies
os.system("cd gatk FindMendelianViolations \
-I unique_variants.recode.vcf.gz --THREAD_COUNT 16 \
          -PED trio_info.ped --MIN_DP 20 --MIN_GQ 20 --MIN_HET_FRACTION 0.4\
          --FEMALE_CHROMS contigs_sex_XZ_scaff.list\
          -O report.mendelian_violation_unique_variants.txt")

# %%
errors = pd.read_csv("report.mendelian_violation_unique_variants.txt",sep="\t",comment="#")
errors['Percentage']=errors['TOTAL_MENDELIAN_VIOLATIONS']/errors['NUM_VARIANT_SITES']

# %%
errors["Mendel_error(%)"]=(errors["TOTAL_MENDELIAN_VIOLATIONS"]/errors["NUM_VARIANT_SITES"]) * 100
errors.describe()

# %%
errors.sum()

# %%
def count_generations_per_family(trios):
    generations = {}
    
    # Assign individuals to initial generation (0)
    current_generation = 1
    for trio in trios:
        for individual in trio:
            generations[individual] = 1
    
    # Update generation numbers based on relationships within families
    while True:
        changed = False
        for trio in trios:
            offspring, sire, dam = trio
            if generations[offspring] == current_generation:
                if sire in generations and generations[sire] == 1:
                    generations[sire] = current_generation + 1
                    changed = True
                if dam in generations and generations[dam] == 1:
                    generations[dam] = current_generation + 1
                    changed = True
        if not changed:
            break
        current_generation += 1
    
    # Create dictionary to store generation count per family
    family_generations = {}
    
    # Calculate generations per family
    for trio in trios:
        family = tuple(sorted(trio))
        if family not in family_generations:
            family_generations[family] = generations[max(trio)]  # Selecting maximum generation within the trio
    
    return family_generations

family_generations = count_generations_per_family(trios)
print(family_generations)


# %%
total_snps = 3169 + 2 - 121 # removing results from one unreliable individual (lower coverage) and the liver duplicate
number_of_generations = 46
total_sequence_length = 2*(3234982288) #wilpena
mutation_rate = total_snps / (total_sequence_length * number_of_generations)
print("Mutation rate per base pair per generation: ", mutation_rate)


# %%
sample_info["Sample ID_dup"]=sample_info["Sample ID_dup"].apply(lambda x: "US_"+str(x))
errors = pd.merge(errors,sample_info[['Sample ID_dup','Total Yield','Tissue Source']],how='left',right_on='Sample ID_dup',left_on='OFFSPRING')
errors["Tissue"]=errors["Tissue Source"].replace({'whole blood':1, 'liver':2, 'lung':3, 'RBC pack':4})
errors[['NUM_DIPLOID_DENOVO','Total Yield','Tissue']].corr()

# %%
def calculate_p_value(column1, column2):
    correlation_coef, p_value = pearsonr(column1, column2)
    return p_value

p_values = pd.DataFrame(index=['NUM_DIPLOID_DENOVO','Total Yield','Tissue'], columns=['NUM_DIPLOID_DENOVO','Total Yield','Tissue'])

for col1 in ['NUM_DIPLOID_DENOVO','Total Yield','Tissue']:
    for col2 in ['NUM_DIPLOID_DENOVO','Total Yield','Tissue']:
        p_values.loc[col1, col2] = calculate_p_value(errors[col1], errors[col2])

print(p_values)


# %%
errors[['NUM_DIPLOID_DENOVO','Total Yield','Tissue Source']]

# %%
pheno.index=pheno["Accession"].apply(lambda x: "US_"+str(x))
errors["Age_father"]=[(pheno.loc[errors.loc[i,"OFFSPRING"].split("_l")[0],"DOB"]- pheno.loc[errors.loc[i,"FATHER"],"DOB"]).days /365.25 for i in errors.index]
errors["Age_mother"]=[(pheno.loc[errors.loc[i,"OFFSPRING"].split("_l")[0],"DOB"]- pheno.loc[errors.loc[i,"MOTHER"],"DOB"]).days /365.25 for i in errors.index]
errors["DENOVO_Mutations_peryear"]=(errors["NUM_DIPLOID_DENOVO"]+errors["NUM_HAPLOID_DENOVO"])/((errors["Age_mother"]+errors["Age_father"])/2)
errors

# %%
errors.sum()

# %%
errors[["DENOVO_Mutations_peryear","NUM_DIPLOID_DENOVO","Age_father","Age_mother"]].corr()

# %%
total_snps =  492.477546 -1.590243 - 14.227160
number_of_generations = 46
total_sequence_length = 2*(3234982288) #wilpena
mutation_rate = total_snps / (total_sequence_length * number_of_generations)
print("Mutation rate per base pair per generation per year: ", mutation_rate)

# %%
total_snps 

# %%
0.001/(mutation_rate*2) # estimation for korv based on LTR divergence in Patric Jern's paper

# %%
0.022/(mutation_rate *2) #phacinbeta according to Patric Jern's paper - and our mutation rate

# %%
0.006/(mutation_rate *2) #phacinbetalike according to Patric Jern's paper - and our mutation rate

# %%
(1/1000)/(mutation_rate*2)

# %%
allel.vcf_to_hdf5('all_contigs_filteredv2.recode.vcf.gz', 'all_contigs_filteredv2.h5', fields='*', overwrite=True)

# %%
import h5py
callset = h5py.File('all_contigs_filtered.h5', mode='r')
callset.keys()

# %%
callset["variants"].keys()

# %%
df_variants =pd.DataFrame()
df_variants["CHROM"]=[char.decode('utf-8') for char in callset["variants/CHROM"]]
df_variants["POS"]=callset["variants/POS"]
df_variants["REF"]=[char.decode('utf-8') for char in callset["variants/REF"]]
df_variants["ALT"]=[char[0].decode('utf-8') for char in callset["variants/ALT"]]
df_variants["numalt"]=callset["variants/numalt"]
df_variants["is_snp"]=callset["variants/is_snp"]
df_variants["DP"]=callset["variants/DP"]
#df_variants["FILTER_PASS"]=[char.decode('utf-8') for char in callset["variants/FILTER_PASS"]]
df_variants["QD"]=callset["variants/QD"]
df_variants["QUAL"]=callset["variants/QUAL"]
df_variants["CR"]=gt.count_called(axis=1)/len(callset["samples"])
df_variants.describe()

# %%
df_variants

# %%
gt = allel.GenotypeArray(callset['calldata/GT'])
gt = gt.compress(df_variants["numalt"]==1,axis=0)
print(len(gt))

# %%
len(df_variants2)

# %%
df_variants2=df_variants[df_variants["numalt"]==1].reset_index(drop=True)
df_variants2[["RAF","AAF"]] = gt.count_alleles(max_allele=1).to_frequencies()
df_variants2.describe()

# %%
df_variants2.groupby("is_snp").count()["POS"]

# %%
df_variants2

# %%
df_variants[df_variants["numalt"]>1].groupby("is_snp").count()["POS"]

# %%
#gt_df=pd.DataFrame(gt.to_n_alt(fill=-1),columns=[sample.decode('utf-8') for sample in callset["samples"]])
for col in df_variants2.columns:
    gt_df[col]=list(df_variants2[col])
gt_df

# %%
#get consensus for duplicates - check concordance 
gt_df['dup1']=[gt_df.loc[i,'dup1_liver'] if gt_df.loc[i,'dup1_liver']==gt_df.loc[i,'dup1_lung'] else -1 for i in gt_df.index]
gt_df['dup2']=[gt_df.loc[i,'dup2_nose'] if gt_df.loc[i,'dup2_nose']==gt_df.loc[i,'dup2_lung'] else -1 for i in gt_df.index]
gt_df['dup3']=[gt_df.loc[i,'dup3_liver'] if gt_df.loc[i,'dup3_liver']==gt_df.loc[i,'dup3_spleen'] else -1 for i in gt_df.index]

# %%
#save to hf5
gt_df.to_hdf('filtered_variants_consensus_duplicates.h5', key='df', mode='w')


