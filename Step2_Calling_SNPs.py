# %%
import pandas as pd
import os
import glob
from multiprocessing import Pool


# %%
#prepare reference genome
reference = 'GCA_030178435.1/GCA_030178435.1_ASM3017843v1_genomic.fna'
#os.system('gatk CreateSequenceDictionary -R '+reference)
#os.system('samtools faidx '+reference)

# %%
#call variants
path = 'path to store SNPs'

def call_variants(bam):
    sample = bam.split('/')[-1].split('_sorted.')[0]
    #os.system('samtools index '+bam)
    command='nohup gatk HaplotypeCaller --java-options \
                "-Xmx20g -XX:ConcGCThreads=2" -R '+reference+' -I '+bam+' \
                -O '+path+sample+'.g.vcf.gz -ERC GVCF 2> '+path+sample+'_haplotypecaller.log'
    os.system(command)
with Pool(3) as p:
    names=glob.glob("path to alignments/*_sorted.bam")
    p.map(call_variants,names)

# %%
global path
global reference

def sample_map(list_files,samples):
    """
        Create a file list with all the samples to be merged into the GenomicsDBImport
    """
    file=open(path+'cohort.sample_map','w')
    for i,sample in enumerate(samples):
        file.write(sample+'\t'+list_files[i]+'\n')
    file.close()
    return path+'cohort.sample_map'
        
def genomicsDBimport(chro=None):
    """
        The GenomicsDBImport tool takes in one or more single-sample GVCFs and 
        imports data over at least one genomics interval (this feature is available 
        in v4.0.6.0 and later and stable in v4.0.8.0 and later), and outputs a 
        directory containing a GenomicsDB datastore with combined multi-sample data. 
        GenotypeGVCFs can then read from the created GenomicsDB directly and output 
        the final multi-sample VCF.
    """
    file = path+'cohort.sample_map'
    cohort=pd.read_csv(file,sep='\t',names=['sample','file'])
    n=len(cohort)
    if chro==None:
        command='nohup gatk --java-options "-Xmx200g" \
        GenomicsDBImport --reader-threads 30 --sample-name-map '+file+' --batch-size '+str(n)+' -L '+path+'intervals.list \
        --genomicsdb-workspace-path '+path+'database/all_contigs --verbosity ERROR --tmp-dir '+path+' 2>'+path+'DBimport.log'
    else:
        command='nohup gatk --java-options "-Xmx20g" \
        GenomicsDBImport --reader-threads 10 --sample-name-map '+file+' --batch-size '+str(n)+' \
        --genomicsdb-workspace-path '+path+'database/'+str(chro)+' -L '+str(chro)+ ' --verbosity ERROR --tmp-dir '+path+' 2>'+path+str(chro)+'_DBimport.log'
    os.system(command)
        
def run_genotypeGVCF(chro):
    output_folder = path +str(chro)
    file=path+'database/'+str(chro)
    command='nohup gatk --java-options \
        "-Xmx200g -XX:ConcGCThreads=30" GenotypeGVCFs -R '+reference+' \
        -V gendb://"'+file+'" -O '+output_folder+'_raw.vcf --verbosity ERROR --tmp-dir '+path+' 2>'+output_folder+'_genotype.log'
    os.system(command)

#generate db and subsequent vcf file (all chromosomes together)
files = sorted(glob.glob(path+"*.g.vcf.gz"))
names = [file.split('/')[-1].split('.')[0] for file in files]
file = sample_map(files,names)
genomicsDBimport()
run_genotypeGVCF('all_contigs')

# %%
#Adding filtering annotation
#os.system('bgzip path_to_all_contigs_raw.vcf')
#os.system('tabix /path_to_all_contigs_raw.vcf.gz')
os.system('nohup gatk VariantFiltration \
-R path_to_GCA_030178435.1/GCA_030178435.1_ASM3017843v1_genomic.fna \
-V path_to_all_contigs_raw.vcf.gz -O path_to_all_contigs_filterINFO.vcf.gz \
--filter-name "Failed_RMSMappingQuality" --filter-expression "MQ < 40.0" --filter-name "Failed_FisherStrand" --filter-expression "FS > 60.0" \
--filter-name "Failed_QualByDepth" --filter-expression "QD < 2.0" --filter-name "Failed_StrandOddsRatio" --filter-expression "SOR > 3.0" \
--filter-name "Failed_MappingQualityRankSumTest" --filter-expression "MQRankSum < -12.5" --filter-name "Failed_ReadPosRankSumTest" \
--filter-expression "ReadPosRankSumTest < -8.0" --genotype-filter-name "Failed_DP" --genotype-filter-expression "DP < 10" \
--set-filtered-genotype-to-no-call True 2> filtering.log')

# %%
#filtering results
os.system('vcftools --gzvcf path_to_all_contigs_filterINFO.vcf.gz \
--out path_to_all_contigs_filteredv2 --recode --recode-INFO-all --max-missing 0.9 --remove-filtered-all --mac 1')

# %%
os.system("bgzip path_to_all_contigs_filteredv2.recode.vcf")
os.system("tabix path_to_all_contigs_filteredv2.recode.vcf.gz")


