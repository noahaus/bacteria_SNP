#stat_summary.py
#author: Noah A. Legall
#date created: June 13th 2019
#date finished:
#dates edited:
#purpose: creation of alignment stats after alignment pipeline is complete.

import sys # use to access arguments
import os # use in order to call commands from terminal script is called in
import subprocess as sp

#this program will create the stats files that data will be extracted from.
#if there are 15 samples studied then there will be 15 rows in the outputted csv.

fastq_dir = sys.argv[1].strip()
bam_dir = sys.argv[2].strip()
stats_dir = sys.argv[3].strip()

#STEP 1: fastq stats generation.
#gunzip everything in the temp_folder
os.chdir(fastq_dir)
os.system("echo $(pwd)")
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
os.system('echo "sample_name,R1_size(full),R2_size(full),Q_ave_R1,Q_ave_R2,R1_ave_read_length,R2_ave_read_length" >> read_stats.csv') #header of read_stats.csv
for i in range(len(R1_list)):
    sample_name = R1_list[i].replace("_R1_001.fastq","")
    output_R1 = R1_list[i].replace("_R1_001.fastq",".R1.stats.txt")
    output_R2 = R2_list[i].replace("_R2_001.fastq",".R2.stats.txt")
    generate_fastq_stats = "fastx_quality_stats -i {} -o {}"
    os.system(generate_fastq_stats.format(R1_list[i],output_R1))
    os.system(generate_fastq_stats.format(R2_list[i],output_R2))
    print("{} and {} created.".format(output_R1,output_R2))
    R1_size = sp.check_output("wc -c {}".format(output_R1),shell=True)
    R2_size = sp.check_output("wc -c {}".format(output_R2))
    Q_ave_R1 = sp.check_output("cat {} | sed \"1d\" | awk '{{sum+=$6}} END{{print sum/NR}}'".format(output_R1),shell=True)
    Q_ave_R2 = sp.check_output("cat {} | sed \"1d\" | awk '{{sum+=$6}} END{{print sum/NR}}'".format(output_R2),shell=True)
    R1_ave_read_length = sp.check_output("awk '{{if(NR%4==2) {{count++; bases += length}} }} END{{print bases/count}}' {}".format(R1_list[i]),shell=True)
    R2_ave_read_length = sp.check_output("awk '{{if(NR%4==2) {{count++; bases += length}} }} END{{print bases/count}}' {}".format(R2_list[i]),shell=True)    
    os.system("echo \"{},{},{},{},{},{},{}\" >> read_stats.csv".format(sample_name,R1_size,R2_size,Q_ave_R1,Q_ave_R2,R1_ave_read_length,R2_ave_read_length))

os.system("mv read_stats.csv -t {}".format(stats_dir))
os.chdir(bam_dir)

os.system('ls | grep ".sorted.bam" > bam_call.txt')
bam = open("bam_call.txt", "r")

bam_list = []

for line in bam:
    bam_list.append(line.strip())

os.system('touch bam_stats.csv')
os.system('echo "total_mapped_reads,ave_coverage,unmapped_reads" >> bam_stats.csv')
for i in range(len(bam_list)):
    ave_coverage = sp.check_output("samtools depth {} | awk '{{sum+=$3}} END {{ print sum/NR}}'".format(bam_list[i]),shell=True)
    total_reads = sp.check_output("samtools flagstat {} | awk -F " " 'NR==1 {{print $1}}'".format(bam_list[i]),shell=True)
    mapped_reads = sp.check_output("samtools flagstat {} | awk -F " " 'NR==1 {{print $5}}'".format(bam_list[i]),shell=True)
    unmapped_reads = total_reads - mapped_reads
    os.system("echo \"{},{},{}\" >> bam_stat.csv".format(mapped_reads,ave_coverage,unmapped_reads))
os.system("mv bam_stats.csv -t {}".format(stats_dir))
os.chdir(stats_dir)
os.system("paste -d \",\" read_stats.csv bam_stats.csv > stat_summary.csv")
