#PBS -S /bin/bash
#PBS -q batch
#PBS -N BWA_align_and_Sort_BAM
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=72:00:00
#PBS -l mem=16gb
#PBS -M noahaus@uga.edu
#PBS -m abe

#cd $PBS_O_WORKDIR
echo "Bash version ${BASH_VERSION}..."
#Let's start by creating the structure of the output folder.
#At the end of the analysis
OUT=$(pwd)output_dir
BAM=$(pwd)output_dir/BAM
BASIC=$(pwd)output_dir/BAM/basic
NODUP=$(pwd)output_dir/BAM/nodup
VCF=$(pwd)output_dir/VCF
FILTER=$(pwd)output_dir/VCF/filtered
PILEUP=$(pwd)output_dir/VCF/pileup
RAW=$(pwd)output_dir/VCF/raw
RAXML=$(pwd)output_dir/RAXML

#Variables for scripts in the package.
STEP_1=$(pwd)/bacteria_SNP/pairread2sortBAM.py
STEP_2=$(pwd)/bacteria_SNP/remove_duplicates.py
STEP_3=$(pwd)/bacteria_SNP/vcf2phylip.py

mkdir $OUT $BAM $BASIC $NODUP $VCF $FILTER $PILEUP $RAW $RAXML

module add BWA/0.7.17-foss-2016b
module add SAMtools/1.9-foss-2016b
#Variable for the reference genome.
#Let's index the reference genome.
REF=$1
EMAIL=$2
#REF=/scratch/noahaus/pipeline_script/fastq_data_set-tb_complex/NC_002945v4.fasta

#STEP 1: ALIGN TO REFERENCE GENOME
python pairread2sortBAM.py $REF
mv *.sorted.bam -t $BASIC
cd $BASIC
echo "Step 1 of pipeline complete" | mail -s "STEP 1: ALIGN TO REFERENCE GENOME" $EMAIL

#STEP 2: REMOVE DUPLICATE READS
module add picard/2.16.0-Java-1.8.0_144
python $STEP_2
mv *.nodup.sorted.bam -t $NODUP
cd $NODUP
echo "Step 2 of pipeline complete" | mail -s "STEP 2: REMOVE DUPLICATE READS" $EMAIL

#STEP 3: VARIANT CALLING
module add BCFtools/1.9-foss-2016b
ls | grep "nodup.sorted.bam" > bam_list.txt
bcftools mpileup -Ou -f $REF -b bam_list.txt > temp.pileup.vcf
mv temp.pileup.vcf -t $PILEUP
bcftools call -Ou --ploidy 1 -mv temp.pileup.vcf > temp.raw.vcf
mv temp.raw.vcf -t $RAW
bcftools filter -s LowQual -e '%QUAL<20 || TYPE="indel"' temp.raw.vcf > output.filter.vcf
python $STEP_3 -i output.filter.vcf
mv output.filter.vcf -t $FILTER
echo "Step 3 of pipeline complete" | mail -s "STEP 3: VARIANT CALLING" $EMAIL

#STEP 4: RAxML TREE GENERATION
module add RAxML/8.2.11-foss-2016b-mpi-avx
cd $FILTER
mpirun raxmlHPC-MPI-AVX -s output.filter.min4.phy -n isolates -m GTRGAMMA -N 100 -p 1000
mv *.53_isolates.*  *.53_isolates -t $RAXML
echo "Step 4 of pipeline complete" | mail -s "STEP 4: RAxML TREE GENERATION" $EMAIL