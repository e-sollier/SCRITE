#!/bin/sh

[ ! -e $1/fastq ] && mkdir $1/fastq
while IFS= read -r line
do
[ ! -f $1/fastq/${line}_1.fastq ] || [ ! -f $1/fastq/${line}_2.fastq ] &&  bsub -W 30 "fastq-dump --split-files $line -O $1/fastq"
done < "$1/SRR_Acc_List.txt"
