# bacteria_SNP (created by Noah A. Legall)
UGA GACRC Sapelo2 specific script to call variants amongst samples of bacteria pathogens.
Download/clone the repository in the working directory where your samples are located and run the program as instructed here.

FASTQ samples -> Trimmed FASTQ, BAM files, VCF files, RAxML tree, Alignment Stats, Genome Annotation file.

running with bash:

_bash ./bacteria_SNP/bacteria_pathogen_SNP.sh /path/to/ref.fa email address_

running on sapelo2 cluster:

_qsub -v "reference=/path/to/reference.fa,email=email@uga.edu" ./bacteria_SNP/bacteria_pathogen_SNP.sh_

To check on your output, just move into the output_dir/ directory!

further documentation on output is found in the file 'BacteriaSNP.docx'
