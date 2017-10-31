import sys
import endtrim
import read_FASTQ
import os
import itertools
BASE = sys.argv[1]
GROUP = sys.argv[2]
out_dict={}
out_dict['TGAATTC']='1'
out_dict['ACTCTGT']='2'
out_dict['CTCTGTA']='3'
out_dict['TAGTACT']='4'
out_dict['GGTCTTG']='5'
out_dict['ATAGTAT']='6'
out_dict['TCCTATC']='7'
out_dict['CATTCGT']='8'
TSO_Handle_hash = {}
for TSO in ['1','2','3','4','5','6','7','8','Unmatched']:
        #TSO_out_dir = out_dir+TSO+"/"
        #os.system('mkdir %s' % TSO_out_dir)
        TSO_Handle_hash[TSO] = [open(BASE.strip()+"_R1_001.paired."+GROUP+"."+TSO+".temp.fastq",'w'),open(BASE.strip()+"_R2_001.paired."+GROUP+"."+TSO+".temp.fastq",'w')]
TSO_Count={}
Unmatched_dict={}
for tso in xrange(1,17,1):
        TSO_Count[tso]=0
TSO_Count['Unmatched']=0
for R1_read,R2_read in itertools.izip(read_FASTQ.read_fastq_with_convert(open(BASE.strip()+"_R1_001.paired."+GROUP+".temp.fastq")),read_FASTQ.read_fastq_with_convert(open(BASE.strip()+"_R2_001.paired."+GROUP+".temp.fastq"))):
        TSO_index= R1_read[1][0:7]
        try:
                TSO=out_dict[TSO_index]
                TSO_Count[int(TSO)]+=1
        except:
                TSO='Unmatched'
                TSO_Count[TSO]+=1
                try:
                        Unmatched_dict[TSO_index]+=1
                except:
                        Unmatched_dict[TSO_index]=1
        R1_corrected = (R1_read[0],R1_read[1][14:],R1_read[2][14:],R1_read[3])
        endtrim.print_fastq(R1_corrected,TSO_Handle_hash[TSO][0],33,0)
        endtrim.print_fastq(R2_read,TSO_Handle_hash[TSO][1],33,0)

