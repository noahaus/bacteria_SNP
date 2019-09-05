#call_snp.py
#author: Noah A. Legall
#date created: July 3rd 2019
#date finished: July 3rd 2019
#dates edited:
#purpose: use freebayes-parallel software to call high quality snps.

print("""The most time consuming of these steps is the process of finding what is a SNP and what is not.
First, we look at every shared position from the isolate sample and then
determine if the depth of nucleotides found in that position from the BAM file constitutes proof of a variant.
Then we have a candidate for a consensus SNP! We record these potential SNPs in a file format known as VCF,
which contains information on genome position and also type of variant displayed. However, the quality of SNP calls varies, so by placing rigid criteria for the type of SNPs we
wish to study through caveats and rules (hard filtering), we end up with a VCF file that has SNPs that the researcher feels comfortable with.

After creating the hard filtered VCFs for each isolate, we can intersect all these
VCF files to only keep the SNPs that are common in all isolates. After every file has been intersected, we merge them into one VCF file separated by isolates,
making it easy to create files that can be used to create phylogenetic trees.
 """)

import sys # use to access arguments
import os # use in order to call commands from the terminal script is called in
import re

ref_genome = sys.argv[1].strip()
pileup = sys.argv[2].strip()
#create a file that lists all the bam files in sorted orderself.
#will be used for indexing and also basecalling shortly thereafter.
os.system('ls | grep "nodup.sorted.bam$" > bam_list.txt')

bam_list = []
bam = open('bam_list.txt','r')

for line in bam:
    bam_list.append(line.strip())

print("indexing the bam files (if necessary)")
os.system("rm samp_list.txt")
os.system("touch samp_list.txt")
for i in range(len(bam_list)):
    output = bam_list[i].replace(".nodup.sorted.bam",".addsample.bam")
    sm = bam_list[i].replace(".nodup.sorted.bam","")
    print("INDEXING: {}".format(bam_list[i]))
    os.system("samtools index {}".format(bam_list[i]))
    print("ADDING READ GROUP: {}".format(bam_list[i]))
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
    print("\n\n\n")
    raw_output = samp_list[i].replace(".addsample.bam",".unfiltered.vcf")
    raw_updated = samp_list[i].replace(".addsample.bam",".mq_fix.vcf")
    print("SNP CALLING: freebayes-parallel chrom_ranges.txt 8 -E -1 -e 1 -u --strict-vcf -f {}  {} > {}".format(ref_genome.strip(),samp_list[i],raw_output))
    os.system("freebayes-parallel chrom_ranges.txt 8 -E -1 -e 1 -u --strict-vcf -f {}  {} > {}".format(ref_genome.strip(),samp_list[i],raw_output))
    #Code inspired by the
    write_fix = open(raw_updated, 'w+')
    with open(raw_output, 'r') as unfiltered:
        for line in unfiltered:
            line = line.strip()
            new_line = re.sub(r';MQM=', r';MQ=', line)
            new_line = re.sub(r'ID=MQM,', r'ID=MQ,', new_line)
            write_fix.write(new_line+"\n")
        write_fix.close()
    filtered_output = samp_list[i].replace(".addsample.bam",".filtered.vcf")
    os.system(r'vcffilter -f "QUAL > 20" %s > %s' % (raw_updated, filtered_output))
    os.system("bgzip {}".format(filtered_output))
    os.system("bcftools index {}.gz".format(filtered_output))
    intersect_isolates = intersect_isolates+" {}.gz".format(filtered_output)
    print("\n\n\n")
#intersect all the VCFs to get common SNPs between all isolates. merge them all into one vcf file.
os.system("bcftools isec -p +{} -n={} {}".format(pileup,len(samp_list),intersect_isolates))
print("INTESECTING ALL VCFs: bcftools isec -p {} -n={} {}".format(pileup,len(samp_list),intersect_isolates))
os.chdir(pileup)
os.system('rm output.merge.vcf.gz output.merge.vcf vcf_list.txt')
os.system('ls | grep ".vcf$" > vcf_list.txt')

vcf_list = []
vcf = open('vcf_list.txt','r')

for line in vcf:
    vcf_list.append(line.strip())

merge_vcf = ""
for i in range(len(vcf_list)):
    os.system("bgzip {}".format(vcf_list[i]))
    os.system("bcftools index {}.gz".format(vcf_list[i]))
    merge_vcf = merge_vcf + " {}.gz".format(vcf_list[i])

print("MERGING THE VCF: bcftools merge {} > output.merge.vcf".format(merge_vcf))
os.system("bcftools merge {} > output.merge.vcf".format(merge_vcf))
