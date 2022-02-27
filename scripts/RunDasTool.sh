#!/bin/bash                              
#$ -l mem_free=20G,h_stack=128m,h_vmem=20G       
#$ -pe smp 4

echo "**** Job starts ****"              
date                                     

echo "**** info ****"                    
echo "User: ${USER}"                     
echo "Job id: ${JOB_ID}"                 
echo "Job name: ${JOB_NAME}"             
echo "Hostname: ${HOSTNAME}"   
echo "Task id: ${SGE_TASK_ID}" 

export PATH="$PATH:/mnt/bin/hmmer/3.1b2/binaries"

# Run for Mats

# Activate conda
conda activate concoct_env

mkdir Mats

# MetaBat

# Generates the MetaBat table needed for DAS Tool
for i in $(ls ../AssemblyMats/MetaBAT/binsMetaBAT.*.fa); do
	# "<<< $i" passes the variable $i as an argument for sed instead of a file 
	BinN=$(sed 's/\.\.\/AssemblyMats\/MetaBAT\/binsMetaBAT\.//g' <<< $i | sed 's/\.fa//g')
	echo "Sample $i and $BinN"
	grep ">" $i | sed 's/>//g' > BinsMatsMetaBat.tsv
	awk -v OFS='\t' -v count="$BinN" '{print $1, "Bin.", count}' BinsMatsMetaBat.tsv >> Mats/MatsMetaBatContigList.tsv
	sed -i 's/Bin.\t/Bin./g' Mats/MatsMetaBatContigList.tsv
done
rm BinsMatsMetaBat.tsv


# CONCOCT

# Generates the CONCOCT table needed for DAS Tool
for i in $(ls ../CONCOCT/Mats/fasta_bins/*.fa); do
	BinN=$(sed 's/\.\.\/CONCOCT\/Mats\/fasta_bins\///g' <<< $i | sed 's/\.fa//g')
	echo "Sample $i and $BinN"
	grep ">" $i | sed 's/>//g' >BinsMatsConcoct.tsv
	awk -v OFS='\t' -v count="$BinN" '{print $1, "Bin.", count}' BinsMatsConcoct.tsv >> Mats/MatsConcoctContigList.tsv
	sed -i 's/Bin.\t/Bin./g' Mats/MatsConcoctContigList.tsv
done
rm BinsMatsConcoct.tsv

# DAS_Tool
cd Mats

ln -s ../../AssemblyMats/final.contigs.2000.fa Mats.fa
DAS_Tool -i MatsConcoctContigList.tsv,MatsMetaBatContigList.tsv -c Mats.fa -l concoct,metabat -o Mats --threads 4 --write_bins 1 

cd ..


echo "**** Job ends ****"
date
