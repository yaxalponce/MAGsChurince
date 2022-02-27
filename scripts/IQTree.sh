#!/bin/bash
#$ -l mem_free=10G,h_stack=128m,h_vmem=10G 
#$ -pe smp 6

echo "**** Job starts ****"
date 

echo "**** info ****"
echo "User: ${USER}" 
echo "Job id: ${JOB_ID}" 
echo "Job name: ${JOB_NAME}" 
echo "Hostname: ${HOSTNAME}" 
echo "Task id: ${SGE_TASK_ID}" 


cd /mnt/data/yaxal/CC/Phylosift

# Run IQ-Tree
iqtree -s 2021_07_31_Bacteria_alignment_gaps.phy --mem 4G -T 6 -B 1000 -bnni -m LG+F+R10

echo "**** Job ends ****"
date
