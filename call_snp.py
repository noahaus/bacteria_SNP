#call_snp.py
#author: Noah A. Legall
#date created: July 3rd 2019
#date finished: July 3rd 2019
#dates edited:
#purpose: use freebayes-parallel software to call high quality snps.

import sys # use to access arguments
import os # use in order to call commands from the terminal script is called in

ref_genome = sys.argv[1].strip()
pileup = sys.argv[2].strip()
#create a file that lists all the bam files in sorted orderself.
#will be used for indexing and also basecalling shortly thereafter.
os.system('ls | grep "nodup.sorted.bam" > bam_list.txt')

bam_list = []
bam = open('bam_list.txt','r')

for line in bam:
    bam_list.append(line.strip())

print("indexing the bam files (if necessary)")
for i in range(len(bam_list)):
    output = bam_list[i].replace(".nodup.sorted.bam",".addsample.bam")
    sm = bam_list[i].replace(".nodup.sorted.bam","")
    print("INDEXING: {}".format(bam_list[i]))
    os.system("samtools index {}".format(bam_list[i]))
    os.system("java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups INPUT={} OUTPUT={} RGLB=lib1 RGPL=illumina RGPU=unit1  RGSM={} USE_JDK_DEFLATER=true USE_JDK_INFLATER=true".format(bam_list[i],output,sm))
    os.system("echo {} >> samp_list.txt".format(output))
    os.system("samtools index {}".format(output))


samp_list = []
samp = open('samp_list.txt','r')

for line in samp:
    samp_list.append(line.strip())

#base calling with freebayes-parallel. hard filtering, keep only SNPs that passed.
intersect_isolates = ""
for i in range(len(samp_list)):
    raw_output = samp_list[i].replace(".addsample.bam",".unfiltered.vcf")
    os.system("freebayes-parallel chrom_ranges.txt 8 -E -1 -e 1 -u --strict-vcf -f {}  {} > {}".format(ref_genome.strip(),samp_list[i],raw_output))
    filtered_output = samp_list[i].replace(".addsample.bam",".filtered.vcf")
    os.system("bcftools filter -s LowQual -i '%QUAL>150 && TYPE=\"snp\" && INFO/DP > 10' {} > {}".format(raw_output,filtered_output))
    pass_output = samp_list[i].replace(".addsample.bam",".only_pass.vcf")
    os.system("grep -v \"LowQual\" {} > {}".format(filtered_output,pass_output))
    os.system("bgzip {}".format(pass_output))
    os.system("bcftools index {}.gz".format(pass_output))
    intersect_isolates = intersect_isolates+" {}.gz".format(pass_output)

#intersect all the VCFs to get common SNPs between all isolates. merge them all into one vcf file.
os.system("bcftools isec -p {} -n={} {}".format(pileup,len(samp_list),intersect_isolates))
os.chdir(pileup)

os.system('ls | grep ".vcf" > vcf_list.txt')

vcf_list = []
vcf = open('vcf_list.txt','r')

for line in vcf:
    vcf_list.append(line.strip())

merge_vcf = ""
for i in range(len(vcf_list)):
    os.system("bgzip {}".format(vcf_list[i]))
    os.system("bcftools index {}.gz".format(vcf_list[i]))
    merge_vcf = merge_vcf + " {}.gz".format(vcf_list[i])

os.system("bcftools merge {} > output.merge.vcf".format(merge_vcf))
