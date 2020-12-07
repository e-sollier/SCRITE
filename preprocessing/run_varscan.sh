#!/bin/sh

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X M; do
	bsub -n 3 -W 18:00 -R "rusage[mem=2000]" "samtools mpileup -d 1000000 -f ../refdata-gex-GRCh38-2020-A/fasta/genome.fa BAM_merged/chr${chr}.bam | java -jar VarScan.v2.4.4.jar mpileup2snp --output-vcf 1 --min-coverage 500 --min-reads2 40  --min-avg-qual 15 --min-var-freq 0.004 --p-value 1 | python filter_variants.py --minREF 300 --minALT 30 > VCF/chr${chr}.vcf"
done
