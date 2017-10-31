import sys
import os
import glob
import numpy
import itertools
folders = sys.argv[1::]
Gene_IDs = []
Cell_IDs = []
counts = []
count_sum = 0
genes_collected = False
for folder in folders:
	count_files = glob.glob(folder+"/*/*/ReadsPerGene.out.tab")
	for count_file in count_files:
		sys.stderr.write(count_file+"\n")
		count_file_ID = count_file.split("/")[-3:-1]
		FASTQ_ID = count_file_ID[0].split("_")
		Cell_ID = "_".join([FASTQ_ID[0],FASTQ_ID[1],FASTQ_ID[2],count_file_ID[1]])
		Cell_IDs.append(Cell_ID)
		file_counts = []
		count_handle = open(count_file)
		unmapped = int(count_handle.readline().strip().split()[2])
		muli = int(count_handle.readline().strip().split()[2])
		nofeature = int(count_handle.readline().strip().split()[2])
		ambig = int(count_handle.readline().strip().split()[2])
		for line in count_handle:
			la = line.strip().split()
			if not genes_collected:
				Gene_IDs.append(la[0])
			file_counts.append(int(la[2]))
			count_sum += int(la[2])
		counts.append(numpy.array(file_counts))
		genes_collected = True
counts = numpy.array(counts,dtype=int)
#sys.stderr.write("%d %d\n" %(count_sum,numpy.sum(counts)))
combined_counts_hash = {}
for count_array,Cell_ID in itertools.izip(counts,Cell_IDs):
	if Cell_ID in combined_counts_hash:
		combined_counts_hash[Cell_ID] += numpy.array(count_array)
	else:
		combined_counts_hash[Cell_ID] = numpy.array(count_array)
combined_Cell_IDs = combined_counts_hash.keys()
combined_counts = [combined_counts_hash[key] for key in combined_Cell_IDs]
combined_counts = numpy.array(combined_counts)
sys.stderr.write("%d %d %d\n" %(count_sum,numpy.sum(counts),numpy.sum(combined_counts)))
#print(counts.shape)
print("\t".join(combined_Cell_IDs))
for gene,gene_values in itertools.izip(Gene_IDs,combined_counts.T):
	print(gene+"\t"+"\t".join(map(str,gene_values)))
