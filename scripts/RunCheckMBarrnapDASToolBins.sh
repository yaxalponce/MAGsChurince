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
export PATH="$PATH:/mnt/bin/bedtools/bedtools2.29.2/bin"

cd /mnt/data/yaxal/CC/DASTool/

for DIR in $(ls -d *); do

	cd /mnt/data/yaxal/CC/DASTool/
	mkdir Barrnap
	mkdir Barrnap/${DIR}

	# Running CheckM
	
	echo "Run CheckM for ${DIR}"
	
	cd ${DIR}
	rm -fr CheckM
	echo "Running CheckM for ${DIR}"
	checkm lineage_wf -t 4 -x fa ${DIR}_DASTool_bins/ CheckM

	# Running Barrnap
	
	echo "Run Barrnap for ${DIR}"
	PATHBARRNAP=/mnt/data/yaxal/CC/DASTool/Barrnap/${DIR}
	

	cd ${DIR}_DASTool_bins/
	for BINFA in $(ls *.fa); do
		BIN=$(basename ${BINFA} .fa)
		/mnt/data/yaxal/Programs/barrnap/bin/barrnap ${BINFA} > ${PATHBARRNAP}/${BIN}_hits.gff3
		/mnt/data/yaxal/Programs/barrnap/bin/barrnap --kingdom arc --lencutoff 0.2 --reject 0.3 --evalue 1e-05 ${BINFA} > ${PATHBARRNAP}/${BIN}_archaea.gff3
		/mnt/data/yaxal/Programs/barrnap/bin/barrnap --kingdom bac --lencutoff 0.2 --reject 0.3 --evalue 1e-05 ${BINFA} > ${PATHBARRNAP}/${BIN}_bacteria.gff3
		
		# Extract fasta sequences from Barrnap 
		# Bins Bacteria and Archaea
		grep "16S" ${PATHBARRNAP}/${BIN}_hits.gff3 | cut -f1,4,5,7 | awk -v OFS='\t'-v bin="$BIN" '{print $1, $2, $3, bin, "1", $4}' > ${PATHBARRNAP}/${BIN}_coordinates.tsv
		bedtools getfasta -fi ${BINFA} -bed ${PATHBARRNAP}/${BIN}_coordinates.tsv -fo ${PATHBARRNAP}/${BIN}_16S.fa.tmp -s -name
		sed -i 's/::.*)//g' ${PATHBARRNAP}/${BIN}_16S.fa.tmp
		awk '/^>/{$0=$0"_"(++i)}1'  ${PATHBARRNAP}/${BIN}_16S.fa.tmp > ${PATHBARRNAP}/${BIN}_16S.fa

		# Bacterial Bins
		grep "16S" ${PATHBARRNAP}/${BIN}_bacteria.gff3 | cut -f1,4,5,7 | awk -v OFS='\t' -v bin="$BIN" '{print $1, $2, $3, bin, "1", $4}' > ${PATHBARRNAP}/${BIN}_bacteria_coordinates.tsv
		bedtools getfasta -fi ${BINFA} -bed ${PATHBARRNAP}/${BIN}_bacteria_coordinates.tsv -fo ${PATHBARRNAP}/${BIN}_bacteria_16S.fa.tmp -s -name
		sed -i 's/::.*)//g' ${PATHBARRNAP}/${BIN}_bacteria_16S.fa.tmp
		awk '/^>/{$0=$0"_"(++i)}1' ${PATHBARRNAP}/${BIN}_bacteria_16S.fa.tmp > ${PATHBARRNAP}/${BIN}_bacteria_16S.fa
		
		# Archaea Bins
		grep "16S" ${PATHBARRNAP}/${BIN}_archaea.gff3 | cut -f1,4,5,7 | awk -v OFS='\t' -v bin="$BIN" '{print $1, $2, $3, bin, "1", $4}' > ${PATHBARRNAP}/${BIN}_archaea_coordinates.tsv
		bedtools getfasta -fi ${PATHBARRNAP}/${BINFA} -bed ${PATHBARRNAP}/${BIN}_archaea_coordinates.tsv -fo ${PATHBARRNAP}/${BIN}_archaea_16S.fa.tmp -s -name
		sed -i 's/::.*)//g' ${PATHBARRNAP}/${BIN}_archaea_16S.fa.tmp
		awk '/^>/{$0=$0"_"(++i)}1' ${PATHBARRNAP}/${BIN}_archaea_16S.fa.tmp > ${PATHBARRNAP}/${BIN}_archaea_16S.fa
		
		rm *.tmp
		
		cd ..
	
	done

done

echo "**** Job ends ****"
date
