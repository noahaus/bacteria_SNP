#PBS -S /bin/bash
#PBS -q batch
#PBS -N BWA_align_and_Sort_BAM
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=72:00:00
#PBS -l mem=16gb
#PBS -M noahaus@uga.edu
#PBS -m abe

#are we working in bash, or the Sapelo2 (University of Georgia specific) cluster?
if [ -z "$PBS_O_WORKDIR" ]
then
  REF=$1 #Variable for the reference genome. first argument
  EMAIL=$2 #Variable for the email you wish to get notifications from. second arguement
  #if you source activate the associated conda environment, then everything should run just fine.
  module add BWA/0.7.17-foss-2016b
  module add SAMtools/1.9-foss-2016b
  module add picard/2.16.0-Java-1.8.0_144
  module add BCFtools/1.9-foss-2016b
  module add freebayes/1.2.0
  module add RAxML/8.2.11-foss-2016b-mpi-avx
  module add FASTX-Toolkit/0.0.14-foss-2016b
  module add Trimmomatic/0.36-Java-1.8.0_144
else
  REF=${reference}  #Variable for the reference genome. first argument
  EMAIL=${email}  #Variable for the email you wish to get notifications from. second arguement
  #add modules the Sapelo2 way
  module add BWA/0.7.17-foss-2016b
  module add SAMtools/1.9-foss-2016b
  module add picard/2.16.0-Java-1.8.0_144
  module add BCFtools/1.9-foss-2016b
  module add freebayes/1.2.0
  module add RAxML/8.2.11-foss-2016b-mpi-avx
  module add FASTX-Toolkit/0.0.14-foss-2016b
  module add Trimmomatic/0.36-Java-1.8.0_144
fi

echo "Bash version ${BASH_VERSION}..."

#Let's start by creating the structure of the output folder.
#At the end of the analysis, all of these directories will have files associated with Variant Calling
#Check out the documentation for the code and the README to see what is in each directory.
OUT=$(pwd)/output_dir
BAM=$(pwd)/output_dir/BAM
BASIC=$(pwd)/output_dir/BAM/basic
NODUP=$(pwd)/output_dir/BAM/nodup
VCF=$(pwd)/output_dir/VCF
FILTER=$(pwd)/output_dir/VCF/filtered
PILEUP=$(pwd)/output_dir/VCF/pileup
RAW=$(pwd)/output_dir/VCF/raw
RAXML=$(pwd)/output_dir/RAXML
FASTQ=$(pwd)/output_dir/FASTQ
STATS=$(pwd)/output_dir/STATS

#Variables for scripts in the package.
STEP_1=$(pwd)/bacteria_SNP/pairread2sortBAM.py
STEP_2=$(pwd)/bacteria_SNP/remove_duplicates.py
STEP_3=$(pwd)/bacteria_SNP/vcf2phylip.py
STEP_5=$(pwd)/bacteria_SNP/stat_summary.py
STEP_6=$(pwd)/bacteria_SNP/snp_table.py
CALL_SNP=$(pwd)/bacteria_SNP/call_snp.py
CHROM=$(pwd)/bacteria_SNP/create_chrom.py

#create the output structure.
mkdir $OUT $BAM $BASIC $NODUP $VCF $FILTER $PILEUP $RAW $RAXML $FASTQ $STATS

#move our fastq samples from the current working directory into the FASTQ folder. this way we always have a way to access our original data
cp *.fastq* -t $FASTQ

#STEP 1: ALIGN TO REFERENCE GENOME
#We should have files separated by R1 and R2. exploit this fact in order to proceed
#with mapping to the reference genome you provide.
cd $FASTQ
python $STEP_1 $REF
mv *.sorted.bam -t $BASIC
echo "Step 1 of pipeline complete" | mail -s "STEP 1: READ QUALITY TRIMMING AND ALIGN TO REFERENCE GENOME" $EMAIL

#STEP 2: REMOVE DUPLICATE READS
#This step ensures that no reads that are the result of PCR duplication and are non informative
#do not complicate further downstream analysis.
cd $BASIC
python $STEP_2
mv *.nodup.sorted.bam -t $NODUP
echo "Step 2 of pipeline complete" | mail -s "STEP 2: REMOVE DUPLICATE READS" $EMAIL

#STEP 3: VARIANT CALLING
#For each sample, use freebayes to call SNPs based on the reference genome.
#perform hard filtering to improve the confidence in SNPs that we find.
cd $NODUP
python $CHROM $REF #We use freebayes-parallel in this pipeline, and in order to do that, we need to break up the reference genome into ranges.
python $CALL_SNP $REF $PILEUP
mv *unfiltered.vcf -t $RAW; mv *.filtered.vcf *only_pass.vcf -t $FILTER;
echo "Step 3 of pipeline complete" | mail -s "STEP 3: VARIANT CALLING" $EMAIL

#STEP 4: RAxML TREE GENERATION
#A basic ML tree for the analysis, can show evolutionary relationship based on nucleotide substitution model.
cd $PILEUP
MERGE=$(pwd)/output.merge.vcf
python $STEP_3 -i $MERGE
PHY=$(ls | grep "merge.*.vcf")
mpirun raxmlHPC-MPI-AVX -s $PHY -n isolates -m GTRGAMMA -N 100 -p 1000
mv *.isolates.*  *.isolates -t $RAXML
echo "Step 4 of pipeline complete" | mail -s "STEP 4: RAxML TREE GENERATION" $EMAIL

#STEP 5: ALIGNMENT STATS
#A file of statistics associated with this run of the pipeline. Let's
#use these stats in order to create summary figures of the reads and genome mapping.
cd $STATS
python $STEP_5 $FASTQ $NODUP $STATS
echo "Step 5 of pipeline complete" | mail -s "STEP 5: ALIGNMENT STATS" $EMAIL

#STEP 6: SNP TABLE (CSV)
#A csv file that shows the isolates on the rows, and the SNP sites (with associated QUAL scores) in the columns.
python $STEP_6 $MALFORM $PASS_ONLY $PHY
echo "Step 6 of pipeline complete" | mail -s "STEP 6: SNP TABLE (CSV)" $EMAIL
