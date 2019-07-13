#pairread2sortBAM.py
#author: Noah A. Legall
#date created: May 9th 2019
#date finished: May 10th 2019
#dates edited: May 21st 2019,
#purpose: automates the process of aligning multiple isolates to a reference genomeself. Will output a directory of sorted BAM files

import sys # use to access arguments
import os # use in order to call commands from terminal script is called in
import re # regular expressions
"""
We start by assuming that you possess reads and also a reliable reference genome to align.
We should align to the genome because otherwise the reads we have are just individual blocks of sequence;
how should they relate to each other? BAM files are used to delineate this relationship between reads
and where the read maps to on the reference genome. After successful mapping, we are not guaranteed that
your reads are in the proper order. For downstream analysis, it is best to sort the mapped reads within the BAM file.
"""

#use this python script in the directory you wish to analyze fastq files in

ref_genome = sys.argv[1].strip()

#index the reference genome
print("INDEXING: "+ ref_genome)
os.system('bwa index '+ref_genome)
print("\n\n\n")
#gunzip everything in the temp_folder
print("UNZIPPING FILES")
os.system('gunzip -r ./')
print("FILES UNZIPPED.")
print("\n\n\n")
#separate the paired reads in the directory
os.system('ls | grep "R1" > R_1.txt')
os.system('ls | grep "R2" > R_2.txt')

#read in the R1.txt & R2.txt to variables R1 and R2 respectively
R1 = open("R_1.txt", "r")
R2 = open("R_2.txt", "r")

R1_list = []
R2_list = []


for R1_line in R1:
    R1_list.append(R1_line.strip())

for R2_line in R2:
    R2_list.append(R2_line.strip())



#in the loop, we can now align each paired read file individually to the reference genome fasta file.
for i in range(len(R1_list)):
    output_pre = re.sub(".R1.*.fastq","",R1_list[i])
    print("aligning "+output_pre.strip()+"files")
    os.system("time java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads 4 -phred33 -trimlog trimlog.txt {} {} {} {} {} {} ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20  MINLEN:25".format(R1_list[i],R2_list[i],output_pre+".R1.trimmed.fastq",output_pre+".unmapped.R1",output_pre+".R2.trimmed.fastq",output_pre+".unmapped.R2"))
    os.system("bwa mem -M -t 4 "+ref_genome+" "+output_pre+".R1.trimmed.fastq "+output_pre+".R2.trimmed.fastq | samtools sort -@4 -o "+output_pre+".sorted.bam > output.quiet.log")

print("BAM FILES CREATED")
os.system("rm R_1.txt R_2.txt output.quiet.log")
print("\n\n\n")
