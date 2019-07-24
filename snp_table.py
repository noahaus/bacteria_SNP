#snp_table.py
#author: Noah A. Legall
#date created: June 24th 2019
#date finished:
#dates edited:
#purpose: create a SNP table with quality scores.

import sys # use to access arguments
import os # use in order to call commands from terminal script is called in

phy = sys.argv[1].strip()
vcf = sys.argv[2].strip()



vcf_p = open(vcf, "r")
vcf_p_list = []

for vcf_lines in vcf_p:
    vcf_p_list.append(vcf_lines)


os.system("touch snp_table.csv")

location = []
location.append(" ")
qual = []
qual.append(" ")

for j in range(len(vcf_p_list)):
    if("\tPASS\t" in vcf_p_list[j]):
        arr = vcf_p_list[j].split()
        print(arr[0]+" "+arr[1])
        location.append(arr[0]+" "+arr[1])
        print(arr[5])
        qual.append(arr[5])

phy_out = open(phy, "r")

phy_list = []
for phy_lines in phy_out:
    phy_list.append(phy_lines)

isolates = []
sequence = []
for k in range(len(phy_list)):
    if(k > 0):
        arr2 = phy_list[k].split()
        isolates.append(arr2[0])
        sequence.append(arr2[1])

line1 = ",".join(map(str, location))
os.system("echo "+line1+" >> snp_table.csv")
line2 = ",".join(map(str, qual))
os.system("echo "+line2+" >> snp_table.csv")

for h in range(len(isolates)):
    line3 = ",".join(map(str, list(sequence[h])))
    record = isolates[h]+","+line3
    os.system("echo "+record+" >> snp_table.csv")
