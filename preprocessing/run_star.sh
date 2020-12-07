#!/bin/sh

input="SRR_Acc_List.txt"
while IFS= read -r line
do
  bsub -n 8 -R "rusage[mem=4000]" "STAR --runThreadN 8 --genomeDir ../refdata-gex-GRCh38-2020-A/star_genome --sjdbGTFfile ../refdata-gex-GRCh38-2020-A/genes/genes.gtf --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix alignments/${line} --readFilesIn fastq/${line}_1.fastq fastq/${line}_2.fastq"
done < "$input"
