# bacteria_SNP (created by Noah A. Legall)
Sapelo2 specific script to call variants amongst samples of bacteria pathogens.
Download/clone the repository in the working directory where your samples are located and run the program as instructed here. 

FASTQ samples -> BAM files, VCF files, RAxML tree, Alignment Stats.

running with bash:

_bash ./bacteria_SNP/bacteria_pathogen_SNP.sh <ref.fa> <email address>_

running on sapelo2 cluster:

_qsub -v "reference=reference.fa,email=email@uga.edu" ./bacteria_SNP/bacteria_pathogen_SNP.sh_

caveats: For pair-read data, please have the suffix of the fastq file be '...R1.fastq'
