#!/bin/bash                              
#$ -l mem_free=10G,h_stack=128m,h_vmem=10G       
#$ -pe smp 1
#$ -t 1:12


echo "**** Job starts ****"              
date                                     

echo "**** info ****"                    
echo "User: ${USER}"                     
echo "Job id: ${JOB_ID}"                 
echo "Job name: ${JOB_NAME}"             
echo "Hostname: ${HOSTNAME}"   
echo "Task id: ${SGE_TASK_ID}" 


cd /mnt/data/yaxal/CC

SAMPLE=$(cat manifest.sediment | awk "NR==${SGE_TASK_ID}")

cd data/Trimmed

# Convert to interleaved sequences

echo "Joining sequences as interleaved"

/mnt/data/yaxal/Programs/bbmap/reformat.sh in1=${SAMPLE}_R1.fq in2=${SAMPLE}_R2.fq out=Interleaved/${SAMPLE}.fq

echo "**** Job ends ****"
date
