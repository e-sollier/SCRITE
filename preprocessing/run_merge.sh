#!/bin/sh

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X M; do
	bsub -W 16:00 -R "rusage[mem=16000]" "samtools merge -r -R chr${chr} -b bam_list.txt BAM_merged2/chr${chr}.bam"
done
