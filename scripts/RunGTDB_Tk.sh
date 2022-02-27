#!/bin/bash                              
#$ -l mem_free=70G,h_stack=128m,h_vmem=70G
#$ -pe smp 7


echo "**** Job starts ****"              
date                                     

echo "**** info ****"                    
echo "User: ${USER}"                     
echo "Job id: ${JOB_ID}"                 
echo "Job name: ${JOB_NAME}"             
echo "Hostname: ${HOSTNAME}"   
echo "Task id: ${SGE_TASK_ID}" 

export PATH="$PATH:/mnt/bin/hmmer/3.1b2/binaries"
export GTDBTK_DATA_PATH=/mnt/data/yaxal/CC/data/GTDB-Tk_release95/

# Run GTDB-Tk 
echo ""
echo "Running GTDB-Tk "

cd cd /mnt/data/yaxalCC/DASTool/

mkdir GTDB_Tk

for DIR in Bordetella Mats Sediment Water; do
	
	cd /mnt/data/yaxal/CC/DASTool/${DIR}
	
	# Running GTDG-Tk
  echo ""
	echo "Running GTDB-Tk for ${DIR}"
        gtdbtk classify_wf --genome_dir ${DIR}_DASTool_bins/ --out_dir /mnt/data/yaxal/CC/DASTool/GTDB_Tk/${DIR} --cpus 7 --extension fa --pplacer_cpus 1 --scratch_dir /mnt/data/yaxal/CC/DASTool/GTDB_Tk/Scratch
 

done

echo "**** Job ends ****"
date
