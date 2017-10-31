import os
import sys
import itertools
import glob
R1 = sys.argv[1]
R2 = sys.argv[2]
R1_trimmed = sys.argv[3]
R2_trimmed = sys.argv[4]
Chunk_Size = 250000
os.system('gunzip -c %s | parallel --pipe -j4 -L4 -N%d \'cat >%s.{#}.temp.fastq\'' % (R1,Chunk_Size,R1) )
os.system('gunzip -c %s | parallel --pipe -j4 -L4 -N%d \'cat >%s.{#}.temp.fastq\'' % (R2,Chunk_Size,R2) )
R1_files = sorted(glob.glob(R1+".*.temp.fastq"), key = lambda x: x.split(".")[-3])
R2_files = sorted(glob.glob(R2+".*.temp.fastq"), key = lambda x: x.split(".")[-3])
if len(R1_files) != len(R2_files):
	exit(0)
os.system('seq %d | parallel --jobs 4 \'cutadapt -m 25 -a %s -A %s -o %s.{}.temp.trimmed.fastq -p %s.{}.temp.trimmed.fastq %s.{}.temp.fastq %s.{}.temp.fastq\'' % (len(R1_files),"file:NexteraPE-PE.fa","file:NexteraPE-PE.fa",R1,R2,R1,R2))

os.system('cat %s.*.temp.trimmed.fastq > %s' % (R1,R1_trimmed))
os.system('cat %s.*.temp.trimmed.fastq > %s' % (R2,R2_trimmed))
os.system('rm %s.*.temp.trimmed.fastq' % (R1))
os.system('rm %s.*.temp.trimmed.fastq' % (R2))
os.system('rm %s.*.temp.fastq' % (R1))
os.system('rm %s.*.temp.fastq' % (R2))
#cat NA12878_chr14.fasta | parallel --block 50M -j16 --recstart '>' --pipe "cat >{#}.fa;blat -maxGap=3 -minMatch=4 -noHead ~/REFERENCES/chr14.fa {#}.fa {#}.psl;rm {#}.fa"
