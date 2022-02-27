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

cd /mnt/data/yaxal/CC/DASTool/

for SAMPLE in Bordetella Mats Sediment Water; do
    
	cd /mnt/data/yaxal/CC/DASTool/${SAMPLE}

	# BWA mapping
	mkdir Mapping
	mkdir coverage
	cd Mapping

	ln -s ../${SAMPLE}_DASTool_bins/*.fa .
	
	for BINFA in $(ls *.fa); do
		BIN=$(basename ${BINFA} .fa)
		
		bwa index ${BINFA}

		echo "Mapping bin ${BIN} to reads with BWA"
		bwa mem -t 14 ${BINFA} -p /mnt/data/yaxal/CC/data/Trimmed/Interleaved/${SAMPLE}.fq.gz > ${BIN}.sam

		echo "Convert SAM to BAM"
		samtools view -bS -F 4 -@14  ${BIN}.sam -o ${BIN}.bam 

		echo "Sort BAM file"
		samtools sort ${BIN}.bam -o ${BIN}.sorted.bam
		rm ${BIN}.bam ${BIN}.sam
 
		# Generate depth file 
		echo "Create depth file for bin ${BIN}"
		jgi_summarize_bam_contig_depths --outputDepth ${BIN}.outDepth.txt --pairedContigs ${BIN}.outPaired.txt ${BIN}.sorted.bam
		
		# Cut columns needed for mmgenome2
		cd ../coverage
		cut -f1,4 ../Mapping/${BIN}.outDepth.txt > ${BIN}_cov
		cd ../Mapping

	done
done

echo "**** Job ends ****"
