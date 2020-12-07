#!/bin/sh

input = "SRR_Acc_List.txt"
while IFS= read -r line
do
  bsub -W 30 "fastq-dump --split-files $line"
done < "$input"
