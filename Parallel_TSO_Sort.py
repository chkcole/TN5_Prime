import sys
import os
import read_FASTQ
import endtrim
import glob
BASE = sys.argv[1]
Chunk_Size = 100000 
R1_paired_trim = BASE.strip()+"_R1_001.paired.fastq"
R2_paired_trim = BASE.strip()+"_R2_001.paired.fastq"
os.system('cat %s | parallel --pipe -L4 -N%d \'cat >%s_R1_001.paired.{#}.temp.fastq\'' % (R1_paired_trim,Chunk_Size,BASE))
os.system('cat %s | parallel --pipe -L4 -N%d \'cat >%s_R2_001.paired.{#}.temp.fastq\'' % (R2_paired_trim,Chunk_Size,BASE))
R1_files = sorted(glob.glob(BASE+"_R1_001.paired.*.temp.fastq"), key = lambda x: x.split(".")[-3])
R2_files = sorted(glob.glob(BASE+"_R2_001.paired.*.temp.fastq"), key = lambda x: x.split(".")[-3])
if len(R1_files) != len(R2_files):
        exit(0)
os.system('seq %d | parallel --jobs 8 \'python TSO_SORT.py %s {#} \'' % (len(R1_files),BASE))
for TSO in ['1','2','3','4','5','6','7','8','Unmatched']:
	os.system('cat '+BASE+"_R1_001.paired.*."+TSO+".temp.fastq > %s"% (BASE.strip()+"_R1_001.paired."+TSO+".fastq"))
	os.system('cat '+BASE+"_R2_001.paired.*."+TSO+".temp.fastq > %s"% (BASE.strip()+"_R2_001.paired."+TSO+".fastq"))
os.system('rm '+BASE+"_R1_001.paired.*.temp.fastq")
os.system('rm '+BASE+"_R2_001.paired.*.temp.fastq")

