#remove_duplicates.py
#author: Noah A. Legall
#date created: May 15th 2019
#date finished: May 15th 2019
#purpose: automates the process of removing duplicates from bam files in a directory.

import sys # use to access arguments
import os # use in order to call commands from terminal script is called in

"""
Duplicate reads are inevitable in the PCR process of creating the sequences we study
(http://www.cureffi.org/2012/12/11/how-pcr-duplicates-arise-in-next-generation-sequencing/
) and if they stay within our BAM file, they could skew later downstream analysis results
by possibly inflating the importance of certain ‘called SNPs’ especially if these variants are
sequencing errors that were replicated. It is best to just remove them all together.
"""
#read in the R1.txt & R2.txt to variables R1 and R2 respectively
os.system('ls | grep ".sorted.bam" > bam_call.txt')
bam = open("bam_call.txt", "r")
os.system('rm bam_call.txt')
bam_list = []

for line in bam:
    bam_list.append(line.strip())

#remove the duplicates one by one
for i in range(len(bam_list)):
    input = bam_list[i]
    output = bam_list[i].replace(".sorted.bam",".nodup.sorted.bam")
    output_pre = bam_list[i].replace(".nodup.sorted.bam","")
    print(input)
    os.system("java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT="+input+" OUTPUT="+output+" METRICS_FILE="+output_pre+".metric.txt REMOVE_SEQUENCING_DUPLICATES=true REMOVE_DUPLICATES=true")


print("non duplicate BAM files created and ready to be moved")
print("DONE")
