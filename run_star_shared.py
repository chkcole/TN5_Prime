import os
import sys
import glob
Data_Dir = sys.argv[1]
fastq_files = glob.glob(Data_Dir+"/*fastq*")
total_jobs = 3
Base_R1 = set()
Base_R2 = set()
unpaired = set()
for fastq_file in fastq_files:
	if "_R1_001" in fastq_file:
		Base_R1.add(fastq_file.split("_R1_001")[0])
	elif "_R2_001" in fastq_file:
		Base_R2.add(fastq_file.split("_R2_001")[0])
	else:
		unpaired.add(fastq_file.split(".fastq")[0])
if Base_R1 != Base_R2:
	print("Unpaired: "+",".join([x for x in Base_R1.union(Base_R2) if x not in Base_R1.intersection(Base_R2)]))

Paired_bases = Base_R1.intersection(Base_R2)
base_list = [[] for i in range(total_jobs)]
for i,base in enumerate(Paired_bases):
	os.system('qsub -q 128-24i.q -cwd -V -v BASE=%s -l mem_free=8G -pe mpi 8 python_TSO_Trimmomatic.sh' % base)
	base_list[i % total_jobs].append(base)
os.system('')
for i,base in enumerate(unpaired):
	os.system('qsub -q 128-24i.q -cwd -V -v BASE=%s -l mem_free=8G -pe mpi 8 python_single_end_no_TSO.sh' % base)
#for i,bases in enumerate(base_list):
#	current_file = "base_list_"+str(i)+".txt"
#	current_handle = open(current_file,'w')
#	for base in bases:
#		current_handle.write(base+"\n")
#	current_handle.close()
#	os.system('qsub -cwd -V -v BASE_LIST=%s -l mem_free=32G -pe mpi 16 Python_STAR_Node.sh' % current_file)
