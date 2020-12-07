import argparse

#split the VCF of the 1000 Genomes Project into different files, one for each chromosome.

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str,help='Input VCF file, containing the data for all chromosomes.')
parser.add_argument('-o', type=str,help='Output directory.')
args = parser.parse_args()

current_chr = "1"
file_chr = open(args.o+ "/chr1.vcf","w")

with open(args.i,"r") as SNP_file:
    for line in SNP_file:
        if line[0]!="#":
            linesplit = line.split("\t")
            chr = linesplit[0]
            if chr!=current_chr:
                file_chr.close()
                print(chr)
                current_chr=chr
                file_chr = open(args.o + "/chr" + str(chr) + ".vcf","w")
            tmp = file_chr.write(line)
file_chr.close()