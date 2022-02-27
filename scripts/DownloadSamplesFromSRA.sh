#!/bin/bash

# Usage bash DownloadSamplesFromSRA.sh FileToRead

while read -r SRA SAMPLE || [[ -n "$sra" ]]; 
do
        printf 'SRA ID: %s, and sample: %s \n' "${SRA}" "${SAMPLE}"
        /mnt/data/yaxal/Programs/sratoolkit.2.9.6-centos_linux64/bin/prefetch -O . ${SRA}
        /mnt/data/yaxal/Programs/sratoolkit.2.9.6-centos_linux64/bin/fastq-dump --split-files --defline-seq '@$sn[_$rn]
/$ri' $SRA.sra
        gzip -c ${SRA}_1.fastq > ${SAMPLE}_R1.fastq.gz 
        gzip -c ${SRA}_2.fastq > ${SAMPLE}_R2.fastq.gz 
        rm -fr ${SRA}*
done < "$1"

