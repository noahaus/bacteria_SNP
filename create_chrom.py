#Quick Script to create a file to use freebayes-parallel functionality.

import sys # use to access arguments
import os # use in order to call commands from terminal script is called in
from Bio import SeqIO

#use this python script in the directory you wish to analyze fastq files in

ref_genome = sys.argv[1].strip()

chrom_ranges = open("chrom_ranges.txt", 'w+')
for record in SeqIO.parse(ref_genome, "fasta"):
    chrom = record.id
    total_len = len(record.seq)
    min_number = 0
    step = 100000
    if step < total_len:
        for chunk in range(min_number, total_len, step)[1:]:
            chrom_ranges.write("{}:{}-{}\n".format(chrom, min_number, chunk))
            min_number = chunk
    chrom_ranges.write("{}:{}-{}\n".format(chrom, min_number, total_len))
chrom_ranges.close()
