#!/bin/sh

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X M; do
	bsub -W 12:00 -n 4 -R "rusage[mem=3000]" "cellSNP -s BAM_merged/chr${chr}.bam -o cellSNP/chr${chr} -b RG_names.txt --cellTAG RG --UMItag None --maxFLAG 255 -R VCF/chr${chr}.vcf"
done
