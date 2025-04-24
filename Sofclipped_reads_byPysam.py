import pysam
import glob
from scipy.stats import mode

def is_partially_mapped(cigar):
    # Check for soft clipping in more than 15 bases
    if (cigar[0][0] == 4 and cigar[0][1] > 20) or (cigar[-1][0] == 4 and cigar[-1][1] > 20): # 4 represents 'Softclipped' in pysam
        return True
    return False

def get_softclipped_sequence_coord(read):
    # Check for soft clipping and return the soft-clipped sequence

    if read.cigar:
        # Check for soft clipping at the beginning of the read
        if read.cigar[0][0] == 4 and read.cigar[0][1] > 20:
            length = int(read.cigar[0][1])
            return (0,length,'left')

        # Check for soft clipping at the end of the read
        elif read.cigar[-1][0] == 4 and read.cigar[-1][1] > 20:
            length = read.cigar[-1][1]
            return (-length,len(read.query_sequence),'right')
    return None


def sam_2_fastq(sampaths, fastq_file):
    """
    Proccess readpairs where at least one read mapped to KoRV (output from samtools view -F 12). 
    Unmapped reads and sofclipped reads are stored in a fastq file for further processing.
    """
    annot = {'gi|425029687|dbj|AB721500.1|':{505:"5'LTR",970:"intergenic",6033:"gag-pol",7894:"env",7936:"intergenic",8500:"3'LTR"},#3'LTR 8440
    'AF151794.2':{505:"5'LTR",961:"intergenic",2526:"gag",2641:"intergenic",6024:"pol",78885:"env",7927:"intergenic",8500:"3'LTR"},#3'LTR 8440
    'gi|516464655|gb|KC779547.1|':{559:"5'LTR",1021:"intergenic",6084:"gag-pol",7966:"env",8008:"intergenic",8600:"3'LTR"},#3'LTR 8566
    "PhER":{478:"5'LTR",6441:"intergentic",7589:"env",7552:"intergenic",8030:"3'LTR"}}#PhER REF: MSTS01000062.1:10912078-10920108 of phaCin_unsw_v4.1
    total = 0; processed = 0
    for samfile in sampaths:
        samfile = pysam.AlignmentFile(samfile, "r", check_sq=False)
        for read in samfile.fetch():
            total+=1
            if read.is_unmapped and read.mate_is_unmapped:
                print(read) # This should never happen -> just checking
                continue
            if not read.is_secondary and not read.is_supplementary: #ony keep primary alignments 
                #(since I have differet strains of the same spp, supplementary alignments are expected, but they can be ignored)
                strand = "R" if read.is_reverse else "F"
                read1 = "Read1" if read.is_read1 else "Read2"
                tstart = read.reference_start #start on KoRV or retrovirus
                tend = read.reference_end #end on KoRV or ...
                if tend==None:
                    tend=tstart+read.qlen
                if read.is_unmapped: #save unmapped or softclipped read to fastq file
                    fastq_file.write("@"+read.qname + "_" +strand+"_"+read1+ "_" + samfile.get_reference_name(read.rname)+
                   "_"+str(tstart)+"-"+str(tend)+ "\n"+ read.query_sequence+"\n+\n"+read.qual+"\n")
                    processed+=1
                elif is_partially_mapped(read.cigar):
                    start,end,side = get_softclipped_sequence_coord(read) 
                    try:
                        if side =='right':
                            for key in annot[samfile.get_reference_name(read.rname)].keys():
                                if tend <= key:
                                    gene =  annot[samfile.get_reference_name(read.rname)][key]
                                    break
                        else:
                            for key in annot[samfile.get_reference_name(read.rname)].keys():
                                if tstart <= key:
                                    gene =  annot[samfile.get_reference_name(read.rname)][key]
                                    break
                    except:
                        gene=""

                            

                    fastq_file.write("@"+read.qname + '_softclipped_'+side+"_" +strand+"_"+read1+'_'+ samfile.get_reference_name(read.rname)+
                     "_"+str(tstart)+"-"+str(tend)+"-"+gene+"\n"+read.query_sequence[start:end] +"\n+\n"+read.qual[start:end]+"\n")
                    processed+=1
    print("TOTAL READS PROCESSED: "+str(total)+"\nREADS INCLUDED: "+str(processed))
    return processed
    
processed = []   
samples = []
for sample in sorted(glob.glob("path to processed data/*_korv_mate_unmapped.sam")):
    sample = sample.split("/")[-1].split("_korv")[0]
    if "L1" in sample:
        sample = sample.split("L1_")[1]
    else:
        sample = sample.split("L2_")[1]
    if sample not in samples:
        samples.append(sample) 
        
for sample in samples:
    fastq_file = open("path to processed data/_"+sample+"_ISflanking_reads.fastq","w")
    print(sample)
    processed.append(sam_2_fastq(glob.glob("path to processed data/*"+sample+"_korv_*.sam"), fastq_file))
    fastq_file.close()