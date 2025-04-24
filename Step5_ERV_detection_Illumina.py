# %%
import glob
import pandas as pd
import gzip
from Bio import SeqIO
import numpy as np
import os
from multiprocessing import Pool

# %%
#Thresholds for the IS detection
threads = 6
min_reads = 20 ## amount of reads to provide evidence for KoRV insertion threshold
global pos_threshold_others
global pos_threshold_korv
pos_threshold_korv=500 ## threshold for distance of reads aligned in genome to identify KoRV insertion site (radius like; because of paired-short-reads)
pos_threshold_others= 9000
mapq_filter = 30
soft_len = 20
reference_koala ="GCA_030178435.1/GCA_030178435.1_ASM3017843v1_genomic.fna.masked"
reference =  'Retroviruses/KoRV_phaCinBeta.fa'


# %%
#get sample info 
sample_info = pd.read_excel('metadata.xlsx')
sample_info

# %%
path = "path to data"
sample_info.index = sample_info["Sample Well"]
names = sorted(glob.glob(path+"*/*R1_001.fastq.gz"))

def trim_and_map2KoRV(file):
    sample = str(sample_info.loc[file.split('/')[6].split('_L')[0],'Sample ID_dup'])
    lane = file.split('/')[6].split('_')[-1]
    folder = path+file.split('/')[6]
    os.system("fastp -i "+file+" -I " + file.split('R1_')[0] +"R2_001.fastq.gz --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --dedup --json "+folder+".json --html "+folder+".html --average_qual 20 \
                -l 35 -p -t 1 --cut_front --cut_tail -q 20 --thread 2 -o "+path+lane+"_"+sample+"R1.fastq.gz -O "+path+lane+"_"+sample+"R2.fastq.gz")

    os.system("bwa-mem2 mem -t 6 -a1 "+reference+" "+path+lane+"_"+sample+"R1.fastq.gz "+path+lane+"_"+sample+"R2.fastq.gz \
        > "+path+lane+"_"+sample+"_korv_mapped.sam" )

    #get one read unmapped and the mate mapped
    os.system("samtools view -f 4 -F 8 -h "+path+lane+"_"+sample+"_korv_mapped.sam > "+path+lane+"_"+sample+"_korv_mate_unmapped.sam")

    #get reads mapped to korv (whose mate is unmapped or mapped)
    os.system("samtools view -F 4 -h "+path+lane+"_"+sample+"_korv_mapped.sam > "+path+lane+"_"+sample+"_korv_allmapped.sam")

    os.system("rm "+path+lane+"_"+sample+"_korv_mapped.sam")
    os.system("rm "+path+lane+"_"+sample+"R1.fastq.gz")
    os.system("rm "+path+lane+"_"+sample+"R2.fastq.gz ")

with Pool(8) as p:
    p.map(trim_and_map2KoRV,names)

# %%
stats = sorted(glob.glob(path+"*bam"))
file = open(path+"list_stats.txt","w")
for stat in stats:
    os.system("samtools idxstats "+stat+" > "+stat+".idxstats")
    file.write(stat+".idxstats \n")
file.close()

# %% [markdown]
# Reads were retrieved by pysam (older Python version = 3.6) in Sofclipped_reads_byPysam.py in pysam environment

# %%
os.system("conda run -n pysam_env python Sofclipped_reads_byPysam.py")

# %%
names = sorted(glob.glob(path+"*.fastq"))
def map_2_koala(file):
    sample = str(file.split('/')[-1].split('_IS')[0])
    ## Align softclipped and unmapped reads to the Koala Genome
    os.system("bwa-mem2 mem -t "+str(threads)+" -a1 "+reference_koala+" "+file+" > "+path+sample+"_koala_mapped.sam")
   
with Pool(4) as p:
    p.map(map_2_koala,names)

# %%
#run IS detection on psam env
os.system("conda run -n pysam_env python integration_detection_byPsam.py")


