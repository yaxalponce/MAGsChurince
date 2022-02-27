#!/bin/bash
#$ -l mem_free=32G,h_stack=128m,h_vmem=32G
#$ -pe smp 1

cd /mnt/data/yaxal/CC/data

mkdir Trimmed

for sample in $(ls *1.fastq.gz)
        do base=$(basename $sample "1.fastq.gz")
        java -jar /mnt/bin/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${base}1.fastq.gz ${base}2.fastq.gz Tri
mmed/${base}1.fq /dev/null Trimmed/${base}2.fq /dev/null ILLUMINACLIP:/mnt/bin/trimmomatic/Trimmomatic-0.36/adapters/Tr
uSeq3-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40
done
