# %%
import pandas as pd
import os

# %%
new_file = open("common_wombat/chrmos1-7.fa","w")
file = open("common_wombat/GCA_028626985.1_vu-2k_genomic.fna","r")
line = file.readline()
while ">J" not in line:
    new_file.write(line)
    line = file.readline()
new_file.close()
file.close()

# %%
os.system("makeblastdb -in common_wombat/chrmos1-7.fa -dbtype nucl -out common_wombat/chrmos1-7")

# %%
os.system('blastn -query GCA_030178435.1/GCA_030178435.1_ASM3017843v1_genomic.fna  \
                    -out GCA_030178435.1/contigs_blast_wombatv2.txt -perc_identity 80 -mt_mode 1 -subject_besthit \
                    -db common_wombat/chrmos1-7 \
                    -outfmt "6 qaccver saccver pident length mismatch gapopen qlen qstart qend sstart send sstrand evalue bitscore" \
                    -num_threads 20 -num_alignments 5 -qcov_hsp_perc 90')

# %%
blast = pd.read_csv("GCA_030178435.1/contigs_blast_wombat.txt", sep='\t', names=["qaccver","saccver",
"pident","length","mismatch","gapopen","qlen","qstart","qend","sstart","send","sstrand","evalue","bitscore"])
blast["q_coverage"]=blast["length"]/blast["qlen"]
blast

# %%
blast["saccver"]=blast["saccver"].replace({"CM053637.1":1,"CM053638.1":2,"CM053639.1":3,"CM053640.1":4,"CM053641.1":5,"CM053642.1":6,"CM053643.1":7})
blast[~blast.astype({"saccver":str})["saccver"].str.contains("JAQ")].sort_values("saccver").groupby("saccver")["qaccver"].count()

# %%
#trying to define chromosomes, based on the number of matches per contig, and then using the first match as position to sort
# Group by qaccver and saccver, summing q_coverage
df_blast = blast.groupby(["qaccver", "saccver"])["q_coverage"].sum().reset_index()

# Dictionaries to store data
chro_contigs = {}
chro_positions = {}
chro_qlen = {}
coverages = {}
chromosomes = {
    "JAOEJA010000856.1":"X","JAOEJA010000095.1":"X","JAOEJA010000096.1":"X","JAOEJA010000141.1":"X",
    "JAOEJA010001005.1":"X","JAOEJA010000416.1":"X","JAOEJA010000058.1":"X","JAOEJA010000813.1":"X",
    "JAOEJA010000812.1":"X","JAOEJA010000139.1":"X","JAOEJA010001021.1":"X","JAOEJA010001063.1":"X",
    "JAOEJA010000249.1":"X","JAOEJA010000059.1":"X","JAOEJA010000503.1":"X","JAOEJA010001019.1":"X",
    "JAOEJA010000532.1":"X","JAOEJA010000412.1":"X","JAOEJA010001020.1":"X","JAOEJA010000385.1":"X",
    "JAOEJA010000808.1":"X","JAOEJA010000520.1":"X","JAOEJA010000356.1":"X","JAOEJA010000415.1":"X",
    "JAOEJA010000373.1":"X","JAOEJA010000817.1":"X","JAOEJA010000335.1":"X","JAOEJA010000413.1":"X"
}

# Process each qaccver
for qaccver, group in df_blast.groupby("qaccver"):
    max_coverage = group["q_coverage"].max()
    max_coverage_saccver = group.loc[group["q_coverage"].idxmax()]["saccver"]
    
    if 'JAQ' in str(max_coverage_saccver) or (qaccver in list(chromosomes)): #keep the assignment of X
        continue
    coverages[qaccver] = max_coverage
    chromosomes[qaccver] = max_coverage_saccver
    
    # Initialize dictionary entries if not present
    if max_coverage_saccver not in chro_contigs:
        chro_contigs[max_coverage_saccver] = []
        chro_positions[max_coverage_saccver] = []
        chro_qlen[max_coverage_saccver] = 0
        
    # Collect data
    chro_contigs[max_coverage_saccver].append(qaccver)
    chro_positions[max_coverage_saccver].append(blast[(blast["qaccver"] == qaccver) & (blast["saccver"] == max_coverage_saccver)]["sstart"].min())
    chro_qlen[max_coverage_saccver] += int(blast[blast["qaccver"] == qaccver].iloc[0]["qlen"])
chro_qlen

# %%
for chro in chro_contigs.keys():
    print(str(chro)+' '+str(len(chro_contigs[chro])))
print('X '+str(list(chromosomes.values()).count('X')))

# %%
import json
with open('GCA_030178435.1/contigs_to_chromosomesv2.json', 'w') as file:
    json.dump(chromosomes, file) 
with open('GCA_030178435.1/contigs_per_chromosome_postionv2.json', 'w') as file:
    json.dump(chro_positions, file) 
with open('GCA_030178435.1/contigs_per_chromosomev2.json', 'w') as file:
    json.dump(chro_contigs, file) 

# %%
#checking for Sex chromosomes
import pandas as pd
sex=pd.read_csv("contigs_sex_scaffSex.tsv",sep="\t")
sex[sex['X_Z']]['Length'].sum()


