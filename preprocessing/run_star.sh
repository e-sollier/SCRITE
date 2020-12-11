#!/bin/sh

[ ! -e $1/BAM ] && mkdir $1/BAM
while IFS= read -r line
do
  [ ! -f $1/BAM/${line}.bam ] && bsub -n 8 -R "rusage[mem=4000]" "STAR --runThreadN 8 --genomeDir refdata-gex-GRCh38-2020-A/star_genome --sjdbGTFfile refdata-gex-GRCh38-2020-A/genes/genes.gtf --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $1/alignments/${line} --readFilesIn $1/fastq/${line}_1.fastq $1/fastq/${line}_2.fastq;mv $1/alignments/${line}Aligned.sortedByCoord.out.bam $1/BAM/${line}.bam; samtools index $1/BAM/${line}.bam"
done < "$1/SRR_Acc_List.txt"
