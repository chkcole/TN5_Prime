import os
import sys
import endtrim
import read_FASTQ
"""
The purpose of this program is to take a FASTA file containing IG transcripts and figure out which Isotype they belong to,
then return a FASTA file containing the IG transcripts with the isotype annotation
The program will be organized as follows:
1. align the IG transcripts to a reference containing IgG,IgA,IgM,IgD and IgE transcripts
2. parse the alignment file and decide which Isotype each transcript belongs to
3. print FASTA file containing IG transcripts with isotype annotation in the header
"""
IG_transcripts = sys.argv[1]
annotated_IG_transcripts = sys.argv[2]
#IG_transcripts_w_isoform = sys.argv[2]
REF = "~/REFERENCES/IMGT_Human_Constant_Database/IMGT_Human_Constant_Regions.fasta"
bwa_out = IG_transcripts+".aligned.sam"
os.system('bwa mem %s %s > %s' %(REF,IG_transcripts,bwa_out))
ID_to_Isotype = {}
for line in open(bwa_out):
	if line[0] != '@':
		la = line.strip().split()
		if int(la[1]) != 4:
			Best_hit = la[2].split("|")[1]
			Best_hit_score = int(la[13].split(":")[-1])
			Alt_hit_score = int(la[14].split(":")[-1])
			if Best_hit_score > Alt_hit_score:
				ID_to_Isotype[la[0]] = Best_hit.split('*')[0]
			else:
				ID_to_Isotype[la[0]] = Best_hit.split('*')[0]
IG_handle = open(annotated_IG_transcripts,'w',)
for ID,sequence,comments in read_FASTQ.read_fasta(open(IG_transcripts)):
	if ID in ID_to_Isotype:
		comments +=" "+ID_to_Isotype[ID]
		new_read = (ID,sequence,[60]*len(sequence)," "+comments)
		endtrim.print_fasta(new_read,IG_handle,"",0,0)
	else:
		comments +=" "+"N/A"
		new_read = (ID,sequence,[60]*len(sequence)," "+comments)
		endtrim.print_fasta(new_read,IG_handle,"",0,0)
