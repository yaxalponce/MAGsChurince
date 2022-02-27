#!/bin/bash                              
#$ -l mem_free=20G,h_stack=128m,h_vmem=20G       
#$ -pe smp 14


echo "**** Job starts ****"              
date                                     

echo "**** info ****"                    
echo "User: ${USER}"                     
echo "Job id: ${JOB_ID}"                 
echo "Job name: ${JOB_NAME}"             
echo "Hostname: ${HOSTNAME}"   
echo "Task id: ${SGE_TASK_ID}" 

export PATH="$PATH:/mnt/bin/hmmer/3.1b2/binaries"

cd /mnt/data/yaxal/CC/data/Trimmed/Interleaved

# The sequences are already interleaved, it was done with 
# the BBMapMats.sh script. We will start merging all 
# samples from the site. 


echo "Merging sequences sequences from the same site"

cat S* > AllMats.fq
rm S*.fq
gzip AllMats.fq

# Metagenome assembly with MEGAHIT

echo "Assembly with MEGAHIT"

cd /mnt/data/yaxal/CC/

/mnt/data/yaxal/Programs/MEGAHIT-1.2.9/bin/megahit --12 data/Trimmed/Interleaved/AllMats.fq.gz --k-list 23,43,63,83,103,123 -t 14 -o /mnt/data/yaxal/CC/AssemblyMats
 
# Select contigs with lenght <= 2K

echo "Filtering by size >= 2Kbs"
cd AssemblyMats
perl ../FilterFastaByLength.pl 2000 final.contigs.fa >final.contigs.2000.fa

# BWA mapping

mkdir Mapping
cd Mapping

ln -s ../final.contigs.2000.fa AllMats.2000.fa
bwa index AllMats.2000.fa


echo "Mapping to reads with BWA"
bwa mem -t 14 AllMats.2000.fa -p /mnt/data/yaxal/CC/data/Trimmed/Interleaved/AllMats.fq.gz > AllMats.sam

  echo "Convert SAM to BAM"
samtools view -bS -F 4 -@14  AllMats.sam -o AllMats.bam 

# MetaBAT

# Sort BAM file
samtools sort AllMats.bam -o AllMats.sorted.bam

cd ..
#mkdir MetaBAT
cd MetaBAT

# Generate depth file neede for MetaBAT
jgi_summarize_bam_contig_depths --outputDepth outDepth.txt --pairedContigs outPaired.txt ../Mapping/*sorted.bam

# Run MetaBAT
metabat -i ../Mapping/AllMats.2000.fa -a outDepth.txt -o binsMetaBAT -t 14 --minCVSum 0 --saveCls -d -v --minCV 0.1 -m 2000

# Create list needed for DASTool
awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' *.fa > concat.fasta
grep ">" concat.fasta > concat.tsv
sed -i 's/\.//g' concat.tsv
sed -i 's/>//g' concat.tsv

# Run CheckM

echo "Running CheckM"
cd ..
mkdir binsDir
cd binsDir
ln -s ../MetaBAT/binsMetaBAT.*.fa .
cd ..
checkm lineage_wf -t 14 -x fa binsDir/ CheckM
rm -fr binsDir


echo "**** Job ends ****"
date
