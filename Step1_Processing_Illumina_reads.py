# %%
import pandas as pd
import os
import glob
from multiprocessing import Pool

# %%
"""
PROCESSING ILLUMINA READS AND ALIGNMENT
"""
#get sample info 
path = 'data path'
sample_info = pd.read_excel(path + 'metadata file.xlsx')
sample_info

# %%
#Run FASTQC
names = sorted(glob.glob(path+"*/*.fastq.gz")) #list all samples
def run_fastqc(file):
    if not os.path.exists(file.split('.')[0]+ '_fastqc.html'):
        os.system('fastqc -q -o ' + path + file.split('/')[6] + ' -t 4 '+ file)
with Pool(20) as p:
    p.map(run_fastqc,names)

# %%
#Trim reads with fastp and map with Dragen
reference = 'GCA_030178435.1'
names = sorted(glob.glob(path+"*/*R1_001.fastq.gz"))
sample_info.index = sample_info["Sample Well"]

 '''
 Generate hash table to reference using:
 dragen-os --build-hash-table true --ht-reference GCF_002099425.1/GCF_002099425.1_phaCin_unsw_v4.1_genomic.fna --output-directory GCF_002099425.1/
 dragen-os --ht-reference GCF_002099425.1/ --ht-uncompress true
 '''

def trim_map(file):
    os.system("fastp -i "+file+" -I " + file.split('R1_')[0] +"R2_001.fastq.gz -o "+file.split('R1_')[0]+"R1_001_FILTERED.fastq.gz \
        -O "+file.split('R1_')[0]+"R2_001_FILTERED.fastq.gz --unpaired1 "+path+file.split('/')[6]+"/singleton.fastq.gz \
            --dedup --unpaired2 "+path+file.split('/')[6]+"/singleton.fastq.gz --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --cut_front --cut_tail -q 20 \
            -R 'fastp_report' -t 1 --average_qual 20 -l 35 -p --thread 4 ")
    lane = file.split('/')[6].split('_')[-1]
    sample = str(sample_info.loc[file.split('/')[6].split('_L')[0],'Sample ID'])
    
    file = file.split('R1_')[0]+"R1_001_FILTERED.fastq.gz"
    if not os.path.exists(path + file.split('/')[6]+ "/"+lane+"_"+sample +"_single.mapping_metrics.csv"):
        #first map single reads
        os.system("dragen-os -r "+reference+" \
            --num-threads 6 --output-directory "+ path + file.split('/')[6] +" --output-file-prefix "+ lane+"_"+sample +"_single -1 "+ path + file.split('/')[6] +"/singleton.fastq.gz \
            --RGID "+lane+" --RGSM "+sample+" --ref-load-hash-bin 1 >> "+ path + file.split('/')[6]+"/dragen_single.log" )
        os.system("samtools sort -@ 4 -m 10G -O BAM -o "+ path + file.split('/')[6]+ "/"+lane+"_"+sample +"_single.bam "+ path + file.split('/')[6]+ "/"+ lane+"_"+sample +"_single.sam ")
        os.system("rm "+ path + file.split('/')[6]+ "/"+ lane+"_"+sample +"_single.sam ")
        os.system("rm "+ path + file.split('/')[6]+ "/singleton.fastq.gz")

        # then map paired-end reads
        os.system("dragen-os -r "+reference+" \
            --num-threads 6 --output-directory "+ path + file.split('/')[6] +" --output-file-prefix "+ lane+"_"+sample +" -1 "+ file +" -2 "+ file.split('R1_')[0] +"R2_001_FILTERED.fastq.gz \
            --RGID "+lane+" --RGSM "+sample+" --ref-load-hash-bin 1 >> "+ path + file.split('/')[6]+"/dragen.log" )
        os.system("samtools sort -@ 4 -m 10G -O BAM -o "+ path + file.split('/')[6]+ "/"+ lane+"_"+sample +".bam "+ path + file.split('/')[6]+ "/"+ lane+"_"+sample +".sam ")
        os.system("rm "+ path + file.split('/')[6]+ "/"+ lane+"_"+sample +".sam ")
        os.system("rm "+file)
        os.system("rm "+file.split('R1_')[0] +"R2_001_FILTERED.fastq.gz")

with Pool(4) as p:
    p.map(trim_map,names)

# %%
#merge bam files from different lanes
names = sorted(glob.glob(path+"*/L1*.bam"))
path = path+'alignments/'
sample_info.index=sample_info['Sample Well']

def merge_bam(file):
    if 'single' not in file:
        sample=sample_info.loc[file.split('/')[-2].split('_L1')[0],"Sample ID_dup"]
        os.system('samtools merge -f -@ 4 '+path+str(sample)+'.bam '+file+' '+file.split('.')[0]+'_single.bam '+file.split('_L1')[0]+'_L2/*.bam ')
        os.system("rm "+file.split('_L1')[0]+'_L1/*.bam ')
        os.system("rm "+file.split('_L1')[0]+'_L2/*.bam ')
        os.system("samtools sort -@ 4 -o "+path+str(sample)+"_sorted.bam "+path+str(sample)+".bam")
        os.system("rm "+path+str(sample)+".bam")
with Pool(3) as p:
    p.map(merge_bam,names)


# %%
#check coverage for samples
def coverage(file):
    os.system('samtools stats '+file+' > '+file+'_stats.txt')
with Pool(20) as p:
    p.map(coverage,sorted(glob.glob(path+"*.bam")))


