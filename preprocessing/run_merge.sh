#!/bin/sh

[ -f $1/bam_list.txt ] && rm $1/bam_list.txt
[ ! -e $1/BAM_merged ] && mkdir $1/BAM_merged
while IFS= read -r line
do
 echo "$1/BAM/${line}.bam" >> $1/bam_list.txt
done < "$1/SRR_Acc_List.txt"

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X M; do
	bsub -W 16:00 -R "rusage[mem=16000]" "samtools merge -r -R chr${chr} -b $1/bam_list.txt $1/BAM_merged/chr${chr}.bam; samtools index $1/BAM_merged/chr${chr}.bam"
done
