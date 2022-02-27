#!/bin/bash                              
#$ -l mem_free=5G,h_stack=128m,h_vmem=5G       
#$ -pe smp 1


echo "**** Job starts ****"              
date                                     

echo "**** info ****"                    
echo "User: ${USER}"                     
echo "Job id: ${JOB_ID}"                 
echo "Job name: ${JOB_NAME}"             
echo "Hostname: ${HOSTNAME}"   
echo "Task id: ${SGE_TASK_ID}" 

export PATH="$PATH:/mnt/bin/hmmer/3.1b2/binaries"
export PATH="$PATH:/mnt/bin/bedtools/bedtools2.30.0/bin"

cd /mnt/data/yaxal/CC/

mkdir Barrnap

for DIR in $(ls -d Assembly*); do
cd /mnt/data/yaxal/CC/Barrnap
        SAMPLE=$(echo ${DIR} | sed -e "s/Assembly//")

	mkdir ${SAMPLE}
	cd ${SAMPLE}

	# Symbolic links to Assemblies and Bins
	ln -s ../../${DIR}/final.contigs.2000.fa ${SAMPLE}.fa
	ln -s ../../${DIR}/MetaBAT/binsMetaBAT.*.fa .

	# Run Barrnap for the assembly sequences >2000 
	/mnt/data/yaxal/Programs/barrnap/bin/barrnap ${SAMPLE}.fa > ${SAMPLE}_hits.gff3
	/mnt/data/yaxal/Programs/barrnap/bin/barrnap --kingdom arc --lencutoff 0.2 --reject 0.3 --evalue 1e-05 ${SAMPLE}.fa > ${SAMPLE}_archaea.gff3
	/mnt/data/yaxal/Programs/barrnap/bin/barrnap --kingdom bac --lencutoff 0.2 --reject 0.3 --evalue 1e-05 ${SAMPLE}.fa > ${SAMPLE}_bacteria.gff3 
	
	# Extract fasta sequences from Barrnap results
	# Sequences >2000
	grep "16S" ${SAMPLE}_hits.gff3 | cut -f1,4,5,7 | awk -v OFS='\t' -v sampl="$SAMPLE" '{print $1, $2, $3, sampl, "1", $4}' > ${SAMPLE}_coordinates.tsv
	bedtools getfasta -fi ${SAMPLE}.fa -bed ${SAMPLE}_coordinates.tsv -fo ${SAMPLE}_16S.fa.tmp -s -name
	sed -i 's/::.*)//g' ${SAMPLE}_16S.fa.tmp
	awk '/^>/{$0=$0"_"(++i)}1'  ${SAMPLE}_16S.fa.tmp > ${SAMPLE}_16S.fa

	# Bacteria
	grep "16S" ${SAMPLE}_bacteria.gff3 | cut -f1,4,5,7 | awk -v OFS='\t' -v sampl="$SAMPLE" '{print $1, $2, $3, sampl, "1", $4}' > ${SAMPLE}_bacteria_coordinates.tsv
	bedtools getfasta -fi ${SAMPLE}.fa -bed ${SAMPLE}_bacteria_coordinates.tsv -fo ${SAMPLE}_bacteria_16S.fa.tmp -s -name
	sed -i 's/::.*)//g' ${SAMPLE}_bacteria_16S.fa.tmp
	awk '/^>/{$0=$0"_"(++i)}1' ${SAMPLE}_bacteria_16S.fa.tmp > ${SAMPLE}_bacteria_16S.fa
	
	# Archaea
	grep "16S" ${SAMPLE}_archaea.gff3 | cut -f1,4,5,7 | awk -v OFS='\t' -v sampl="$SAMPLE" '{print $1, $2, $3, sampl, "1", $4}' > ${SAMPLE}_archaea_coordinates.tsv
	bedtools getfasta -fi ${SAMPLE}.fa -bed ${SAMPLE}_archaea_coordinates.tsv -fo ${SAMPLE}_archaea_16S.fa.tmp -s -name
	sed -i 's/::.*)//g' ${SAMPLE}_archaea_16S.fa.tmp
	awk '/^>/{$0=$0"_"(++i)}1' ${SAMPLE}_archaea_16S.fa.tmp > ${SAMPLE}_archaea_16S.fa
	
	rm *.tmp
	
	for BINFA in $(ls binsMetaBAT.*.fa); do
		BIN=$(basename ${BINFA} .fa | sed 's/binsMetaBAT\.//g')
		/mnt/data/yaxal/Programs/barrnap/bin/barrnap ${BINFA} > Bin_${BIN}_hits.gff3
		/mnt/data/yaxal/Programs/barrnap/bin/barrnap --kingdom arc --lencutoff 0.2 --reject 0.3 --evalue 1e-05 ${BINFA} > Bin_${BIN}_archaea.gff3
		/mnt/data/yaxal/Programs/barrnap/bin/barrnap --kingdom bac --lencutoff 0.2 --reject 0.3 --evalue 1e-05 ${BINFA} > Bin_${BIN}_bacteria.gff3
		
		# Extract fasta sequences from Barrnap 
		# Bins Bacteria and Archaea
		grep "16S" Bin_${BIN}_hits.gff3 | cut -f1,4,5,7 | awk -v OFS='\t' -v sampl="$SAMPLE" -v bin="$BIN" '{print $1, $2, $3, sampl, ".Bin_", bin, "1", $4}' | sed 's/\t\.Bin_\t/.Bin_/g' >Bin_${BIN}_coordinates.tsv
		bedtools getfasta -fi ${BINFA} -bed Bin_${BIN}_coordinates.tsv -fo Bin_${BIN}_16S.fa.tmp -s -name
		sed -i 's/::.*)//g' Bin_${BIN}_16S.fa.tmp
		awk '/^>/{$0=$0"_"(++i)}1'  Bin_${BIN}_16S.fa.tmp > Bin_${BIN}_16S.fa

		# Bacterial Bins
		grep "16S" Bin_${BIN}_bacteria.gff3 | cut -f1,4,5,7 | awk -v OFS='\t' -v sampl="$SAMPLE" -v bin="$BIN" '{print $1, $2, $3, sampl, ".Bin_", bin, "1", $4}' | sed 's/\t\.Bin_\t/.Bin_/g' >Bin_${BIN}_bacteria_coordinates.tsv
		bedtools getfasta -fi ${BINFA} -bed Bin_${BIN}_bacteria_coordinates.tsv -fo Bin_${BIN}_bacteria_16S.fa.tmp -s -name
		sed -i 's/::.*)//g' Bin_${BIN}_bacteria_16S.fa.tmp
		awk '/^>/{$0=$0"_"(++i)}1' Bin_${BIN}_bacteria_16S.fa.tmp > Bin_${BIN}_bacteria_16S.fa
		
		# Archaea Bins
		grep "16S" Bin_${BIN}_archaea.gff3 | cut -f1,4,5,7 | awk -v OFS='\t' -v sampl="$SAMPLE" -v bin="$BIN" '{print $1, $2, $3, sampl, ".Bin_", bin, "1", $4}' | sed 's/\t\.Bin_\t/.Bin_/g' >Bin_${BIN}_archaea_coordinates.tsv
		bedtools getfasta -fi ${BINFA} -bed Bin_${BIN}_archaea_coordinates.tsv -fo Bin_${BIN}_archaea_16S.fa.tmp -s -name
		sed -i 's/::.*)//g' Bin_${BIN}_archaea_16S.fa.tmp
		awk '/^>/{$0=$0"_"(++i)}1' Bin_${BIN}_archaea_16S.fa.tmp > Bin_${BIN}_archaea_16S.fa
		
		rm *.tmp
	
	done
	
#	cat Bin_+([[:digit:]])_16S.fa > ${SAMPLE}_Bin.fa
	cat Bin*bacteria_16S.fa > ${SAMPLE}_bacteria_Bin.fa 
	cat Bin*archaea_16S.fa > ${SAMPLE}_archaea_Bin.fa
	
done

echo "**** Job ends ****"
date
