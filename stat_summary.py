#stat_summary.py
#author: Noah A. Legall
#date created: June 13th 2019
#date finished:
#dates edited:
#purpose: creation of alignment stats after alignment pipeline is complete.

import sys # use to access arguments
import os # use in order to call commands from terminal script is called in

#this program will create the stats files that data will be extracted from.
#if there are 15 samples studied then there will be 15 rows in the outputted csv.

#STEP 1: fastq stats generation.
#gunzip everything in the temp_folder
print("unzipping files if needed")
os.system('gunzip -r ./')
print("files unzipped")

#separate the pair reads into two files
os.system('ls | grep "R1" > R_1.txt')
os.system('ls | grep "R2" > R_2.txt')

#let's put the file names in a list data structure.
R1 = open("R_1.txt", "r")
R2 = open("R_2.txt", "r")

R1_list = []
R2_list = []


for R1_line in R1:
    R1_list.append(R1_line.strip())

for R2_line in R2:
    R2_list.append(R2_line.strip())

#let's start with some basic read stats
os.system('touch read_stats.csv') #create the read_stats file
os.system('echo "sample_name,R1_size,R2_size,Q_ave_R1,Q_ave_R2,ave_read_length" >> read_stats.csv') #header of read_stats.csv
for i in range(len(R1_list)):
    sample_name = R1_list[i].replace("_R1_001.fastq","")
    output_R1 = R1_list[i].replace("_R1_001.fastq",".R1.stats.txt")
    output_R2 = R2_list[i].replace("_R2_001.fastq",".R2.stats.txt")
    generate_fastq_stats = "fastx_quality_stats -i {} -o {}"
    os.system(generate_fastq_stats.format(R1_list[i],output_R1))
    os.system(generate_fastq_stats.format(R2_list[i],output_R2))
    print("{} and {} created.".format(output_R1,output_R2))
