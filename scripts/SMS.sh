#!/bin/bash
#$ -l mem_free=25G,h_stack=128m,h_vmem=25G 
#$ -pe smp 4


echo "**** Job starts ****"
date 

echo "**** info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}" 
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}" 


cd /mnt/data/yaxal/CC/Phylosift

# Run SMS
bash /mnt/data/yaxal/Programs/sms-1.8.4/sms.sh -i 2021_07_31_Bacteria_alignment_gaps.phy -d aa -c AIC -s NNI -p SMSbacteria

echo "**** Job ends ****"
date
