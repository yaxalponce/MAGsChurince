#!/bin/bash                              
#$ -l mem_free=10G,h_stack=128m,h_vmem=10G       
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

cd /mnt/data/yaxal/CC/
mkdir CONCOCT

for DIR in $(ls -d Assembly*); do
	conda activate concoct_env
	
	cd /mnt/data/yaxal/CC/CONCOCT
        SAMPLE=$(echo ${DIR} | sed -e "s/Assembly//")

	mkdir ${SAMPLE}
	cd ${SAMPLE}

	ln -s ../../${DIR}/final.contigs.2000.fa ${SAMPLE}.fa
	ln -s ../../${DIR}/MetaBAT/outDepth.txt .

	echo "Running CONCOCT for sample ${SAMPLE}"
	concoct --coverage_file outDepth.txt --composition_file ${SAMPLE}.fa --threads 4 -b ${SAMPLE}_concoct
	merge_cutup_clustering.py ${SAMPLE}_concoct_clustering_gt1000.csv  > ${SAMPLE}_concoct_clustering_merged.csv

	mkdir fasta_bins
	extract_fasta_bins.py ${SAMPLE}.fa ${SAMPLE}_concoct_clustering_merged.csv --output_path fasta_bins

	# Run CheckM
	conda deactivate
	echo "Running CheckM"
	checkm lineage_wf -x fa fasta_bins/ CheckM
	
	cat fasta_bins/*.fa > fasta_bins/Concoct_${SAMPLE}_Bins_cat.fas
	
done
	
	
echo "**** Job ends ****"
date
