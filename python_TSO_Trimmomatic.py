import sys
import os
import itertools
import read_FASTQ
import endtrim
base = sys.argv[1]
out_dict={}
out_dict['TGAATTC']='1'
out_dict['ACTCTGT']='2'
out_dict['CTCTGTA']='3'
out_dict['TAGTACT']='4'
out_dict['GGTCTTG']='5'
out_dict['ATAGTAT']='6'
out_dict['TCCTATC']='7'
out_dict['CATTCGT']='8'
out_dict['ATTATCC']='9'
#out_dict['AACCATT']='10'
#out_dict['TACATCA']='11'
#out_dict['TGTTCTA']='12'
#out_dict['TCGCCAT']='13'
#out_dict['GTGTGAC']='14'
#out_dict['TGACGCT']='15'
#out_dict['TATCGTC']='16'
R1 = base.strip()+"_R1_001.fastq.gz"
R1_paired_trim =base.strip()+"_R1_001.paired.fastq"
R1_unpaired_trim =base.strip()+"_R1_001.unpaired.fastq"
R2 = base.strip()+"_R2_001.fastq.gz"
R2_paired_trim = base.strip()+"_R2_001.paired.fastq"
R2_unpaired_trim =base.strip()+"_R2_001.unpaired.fastq"
out_dir = base.strip()+"/"
trimlog = out_dir+"trimmomatic_log.txt"
os.system('mkdir %s' % out_dir)

#R1_paired_trim_IG_Filtered=base.strip()+"_R1_001.paired.IG_Filtered.fastq"
#R2_paired_trim_IG_Filtered=base.strip()+"_R2_001.paired.IG_Filtered.fastq"

#ADAPTERS
#PrefixNX_1 = "AGATGTGTATAAGAGACAG"
#PrefixNX_2 = "AGATGTGTATAAGAGACAG"
#Trans1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
#Trans1_rc = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"
#Trans2 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
#Trans2_rc = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
#Adapters = [PrefixNX_1,PrefixNX_2,Trans1,Trans1_rc,Trans2,Trans2_rc]

#os.system('java -jar /afs/cats.ucsc.edu/users/q/chkcole/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 %s %s %s %s %s %s ILLUMINACLIP:/afs/cats.ucsc.edu/users/q/chkcole/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > %s' % (R1,R2,R1_paired_trim,R1_unpaired_trim,R2_paired_trim,R2_unpaired_trim,trimlog) )



#os.system("~/.local/bin/cutadapt -m 25 -a %s -A %s -o %s -p %s %s %s" %("file:NexteraPE-PE.fa","file:NexteraPE-PE.fa",R1_paired_trim,R2_paired_trim,R1,R2))
os.system("python Parallel_CutAdapt.py %s %s %s %s" %(R1,R2,R1_paired_trim,R2_paired_trim))

TSO_Handle_hash = {}
for TSO in ['1','2','3','4','5','6','7','8','Unmatched']:
	TSO_out_dir = out_dir+TSO+"/"
	os.system('mkdir %s' % TSO_out_dir)
#	TSO_Handle_hash[TSO] = [open(base.strip()+"_R1_001.paired."+TSO+".fastq",'w'),open(base.strip()+"_R2_001.paired."+TSO+".fastq",'w')]
os.system('python Parallel_TSO_Sort.py %s' %(base.strip()))




#TSO_Count={}
#Unmatched_dict={}
#for tso in xrange(1,17,1):
#	TSO_Count[tso]=0
#TSO_Count['Unmatched']=0
#for R1_read,R2_read in itertools.izip(read_FASTQ.read_fastq_with_convert(open(R1_paired_trim)),read_FASTQ.read_fastq_with_convert(open(R2_paired_trim))):
#	TSO_index= R1_read[1][0:7]
#	try:
#		TSO=out_dict[TSO_index]
#		TSO_Count[int(TSO)]+=1
#	except:
#		TSO='Unmatched'
#		TSO_Count[TSO]+=1
#		try:
#			Unmatched_dict[TSO_index]+=1
#		except:
#			Unmatched_dict[TSO_index]=1
#	R1_corrected = (R1_read[0],R1_read[1][14:],R1_read[2][14:],R1_read[3])
#	endtrim.print_fastq(R1_corrected,TSO_Handle_hash[TSO][0],33,0)
#	endtrim.print_fastq(R2_read,TSO_Handle_hash[TSO][1],33,0)
#os.system('rm %s %s %s %s' %(R1_paired_trim,R1_unpaired_trim,R2_paired_trim,R2_unpaired_trim))




for TSO in ['1','2','3','4','5','6','7','8']:
	#TSO_Handle_hash[TSO][0].close()
	#TSO_Handle_hash[TSO][1].close()
	R1_TSO= base.strip()+"_R1_001.paired."+TSO+".fastq"
	R2_TSO= base.strip()+"_R2_001.paired."+TSO+".fastq"
	R1_TSO_filtered = base.strip()+"_R1_001.paired.IG_Filtered."+TSO+".fastq"
	R2_TSO_filtered = base.strip()+"_R2_001.paired.IG_Filtered."+TSO+".fastq"
	out_dir_TSO= out_dir+TSO+"/"
	#test_of_IG_quantify_folder = out_dir_TSO+"/IG_Quantify/"
	#Transcriptome generation
	os.system("~/SPAdes-3.10.1/bin/rnaspades.py --sc -1 %s -2 %s -o %s" %(R1_TSO,R2_TSO,out_dir_TSO))
	#os.system("~/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --max_memory 10G --CPU 8 --full_cleanup  --left %s --right %s --output %s" %(R1_TSO,R2_TSO, out_dir_TSO+"/trinity_out"))
	#os.system("mv %s %s" % (out_dir_TSO+"/trinity_out.Trinity.fasta", out_dir_TSO+"/transcripts.fasta"))
	#/Transcriptome generation
	#ig Annotation
	os.system("/afs/cats.ucsc.edu/users/q/chkcole/ncbi-igblast-1.6.1/bin/igblastn -germline_db_V ~/REFERENCES/igblast_human_sequences/IMGT_V_formated -germline_db_J ~/REFERENCES/igblast_human_sequences/IMGT_J_formated -germline_db_D ~/REFERENCES/igblast_human_sequences/IMGT_D_formated -domain_system imgt -query %s -out %s -num_alignments_J=1 -num_alignments_D=1 -num_alignments_V=1 -outfmt 7 -evalue=1e-50" % (out_dir_TSO+"/"+"transcripts.fasta",out_dir_TSO+"/"+"igblast_out.txt"))
	os.system("python Extract_IG_Transcripts.py %s %s %s" % (out_dir_TSO+"/"+"igblast_out.txt",out_dir_TSO+"/"+"transcripts.fasta",out_dir_TSO+"/"+"IG_Transcripts.fasta"))
	os.system("python Annotate_IG_W_isoform.py %s %s" % (out_dir_TSO+"/"+"IG_Transcripts.fasta",out_dir_TSO+"/"+"IG_Transcripts_w_isoform.fasta"))
	os.system("python Quantify_isoforms.py %s %s %s %s %s %s" % (R1_TSO,R2_TSO,R1_TSO_filtered,R2_TSO_filtered,out_dir_TSO+"/"+"IG_Transcripts_w_isoform.fasta",out_dir_TSO+"/"+"IG_Table.tab"))
	#/ig Annotation
	#IG_Filtering
	#R1_TSO_filtered = base.strip()+"_R1_001.paired.IG_Filtered."+TSO+".fastq"
	#R2_TSO_filtered = base.strip()+"_R2_001.paired.IG_Filtered."+TSO+".fastq"
	#os.system('bwa mem -t 8 -k 15 -U 0 ~/REFERENCES/BWA_GENCODE_IG_Index/IG_Transcripts.fa %s %s | samtools fastq -F 256 -f 12  -1 %s -2 %s -' % (R1_TSO,R2_TSO,R1_TSO_filtered,R2_TSO_filtered) )
	os.system('~/STAR/bin/Linux_x86_64/STAR --limitBAMsortRAM 10000000000 --outSAMunmapped Within KeepPairs --outSAMtype BAM Unsorted SortedByCoordinate --quantMode GeneCounts --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNoverReadLmax 0.04 --sjdbScore 1 --runThreadN 16 --genomeLoad LoadAndKeep --genomeDir ~/REFERENCES/ensembl_hg38_mRNA --outSAMstrandField intronMotif --readFilesIn %s %s --outFileNamePrefix %s' % (R1_TSO,R2_TSO,out_dir_TSO))
	#os.system('mkdir %s' %(test_of_IG_quantify_folder))
	#os.system('~/STAR/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir ~/REFERENCES/STAR_hg38 --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate  --readFilesIn %s %s --outFileNamePrefix %s' % (R1_TSO,R2_TSO,test_of_IG_quantify_folder))
	os.system('samtools index %s' % (out_dir_TSO+"/Aligned.sortedByCoord.out.bam"))
	os.system('samtools view -L IGH.bed -F 4 %s | wc -l > %s' % (out_dir_TSO+"/Aligned.sortedByCoord.out.bam",out_dir_TSO+"IGH.txt"))#14:105566277-106879844
	os.system('samtools view -L IGK.bed -F 4 %s | wc -l > %s' % (out_dir_TSO+"/Aligned.sortedByCoord.out.bam",out_dir_TSO+"IGK.txt"))#2:88857361-90235368
	os.system('samtools view -L IGL.bed -F 4 %s | wc -l > %s' % (out_dir_TSO+"/Aligned.sortedByCoord.out.bam",out_dir_TSO+"IGL.txt"))#22:22026076-22922913
	#os.system('rm -r %s' %(test_of_IG_quantify_folder))
	#os.system('kallisto quant -t 8 --fr-stranded --index=/afs/cats.ucsc.edu/users/q/chkcole/REFERENCES/KALLISTO_INDEX/gencode.v26.transcripts --output-dir=%s %s %s' % (out_dir_TSO,R1_TSO_filtered,R2_TSO_filtered))
	#os.system('rm %s %s' %(R1_TSO,R2_TSO))
	#os.system('rm %s %s' %(R1_TSO_filtered,R2_TSO_filtered))
	#feature_counts_TSO = out_dir_TSO+"feature_counts.txt"
	#alignments_TSO = out_dir_TSO+"Aligned.out.bam"
	#os.system('~/subread-1.5.1-source/bin/featureCounts -p -R -T 5 -t exon -g gene_id -a ~/REFERENCES/GENCODE/gencode.v26.primary_assembly.annotation.gtf -o %s %s' %(feature_counts_TSO,alignments_TSO))
