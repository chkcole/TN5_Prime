import os
import sys
import glob
import collections
Cell_IDs = set()
Cell_IDs_to_IG = collections.defaultdict(list)
for folder in sys.argv[1::]:
	for IG_Table in glob.glob(folder+"/*/*/IG_Table.tab"):
		#print(IG_Table)
		IG_file_ID = IG_Table.split("/")[-3:-1]
		FASTQ_ID = IG_file_ID[0].split("_")
		Cell_ID = "_".join([FASTQ_ID[0],FASTQ_ID[1],FASTQ_ID[2],IG_file_ID[1]])
		Cell_IDs.add(Cell_ID)
		for line in open(IG_Table):
			la = line.strip().split()
			Cell_IDs_to_IG[Cell_ID].append(la)
for Cell_ID in Cell_IDs:
	for IG_list in Cell_IDs_to_IG[Cell_ID]:
		print(Cell_ID+"\t"+"\t".join(IG_list))

#print(str(len(Cell_IDs))+"\t"+str(len([Cell_ID for Cell_ID in Cell_IDs if len(Cell_IDs_to_IG[Cell_ID]) > 0])))

