import os
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str,help='Input directory containing the outputs of cellSNP (uncompressed).')
parser.add_argument('-o', type=str,help='Output directory.')
parser.add_argument('--freq', type=str,default = "genome1K",help='Directory containing VCF files with the frequency of SNPs in the population.')
parser.add_argument('--minPer', type=float,default = 10.0,help='Minimum Percentage of cells expressing a locus, for this locus to be selected.')
args = parser.parse_args()


def read_1K(file):
    """
    Read VCF file from the 1000 Genomes Project and store the frequency of the variants.
    """
    SNP_freq = {}
    if not os.path.isfile(file):
        return SNP_freq
    with open(file,"r") as SNP_file:
        for line in SNP_file:
            linesplit = line.split("\t")
            AFs = linesplit[7].split(";")[1][3:].split(",")
            alts = linesplit[4].split(",")
            for index_alt in range(len(alts)):
                SNP_freq["chr"+linesplit[0] + "_" + linesplit[1] + "_" + alts[index_alt]] = float(AFs[index_alt])
    return SNP_freq

def keep_locus(ref,alt,oth,minPercentageCellsExpressed=10):
    """
    Returns a boolean indicating whether or not we keep this locus
    """
    reads_alt = np.sum(alt)
    reads_ref = np.sum(ref)
    reads_oth = np.sum(oth)
    if reads_alt >= 6 *reads_oth and reads_ref >=6*reads_oth:
        #filter out loci for which there are many "other" variants. It is position where sequencing errors are common.
        n_cells = len(ref)
        cells_expressed = np.sum((alt+ref)>0)
        cells_alt_expressed = np.sum(alt>0)
        cells_ref_expressed = np.sum(ref>0)
        if cells_expressed / float(n_cells) >= minPercentageCellsExpressed / 100.0:
            # Keep only loci which are expressed in at least some percentage of the cells
            if cells_alt_expressed>= 6 and cells_ref_expressed / float(n_cells)>= 0.7*minPercentageCellsExpressed / 100.0:
                cells_alt_noref = np.sum((alt>4) & (ref==0)) # number of cells where the alt allele is expressed but not the ref allele
                cells_ref_noalt = np.sum((ref>4) & (alt==0)) # number of cells where the ref allele is expressed but not the alt allele
                if cells_alt_noref>=1.2 * cells_ref_noalt or cells_ref_noalt>=1.2 * cells_alt_noref:
                    # Keep only loci for which some cells might be homozygous for one of the alleles.
                    return True
    return False


#Parse output of cellSNP

REFs = []
ALTs = []
loci = []
SNPs =[]

chromosomes_all = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",\
    "chr19","chr20","chr21","chr22","chrX","chrM"] 
chromosomes = []
for x in chromosomes_all:
    if x in os.listdir(args.i):
        chromosomes.append(x)

for filename in chromosomes: 
    print("Reading chromosome " + filename)
    if args.freq is not None:
        SNP_freq = read_1K(os.path.join(args.freq,filename+".vcf"))
    else:
        SNP_freq = {}

    #Parse the output of cellSNP
    with open(os.path.join(args.i,filename),"r") as VCF: 
        line = VCF.readline()
        while line[:2]=="##" or line=="\n":
            line = VCF.readline()
        linesplit = line.split("\t")
        barcodes = linesplit[9:]
        if "SRR" in barcodes[0]:
            barcodes = [barcode[:10] for barcode in barcodes]
        barcodes = [barcode.strip("\n") for barcode in barcodes]

        for line in VCF:
            if line=="\n":
                continue
            linesplit = line.split("\t")
            locus = linesplit[0] + "_" + linesplit[1]
            allele = locus + "_" + linesplit[4]
            
            ref_locus = []
            alt_locus = []
            oth_locus = []
            for cell_info in linesplit[9:]:
                infosplit = cell_info.split(":")
                AD = 0 if infosplit[1]=="." else int(infosplit[1])
                DP = 0 if infosplit[2]=="." else int(infosplit[2])
                OTH = 0 if infosplit[3]=="." else int(infosplit[3])
                alt_locus.append(AD)
                ref_locus.append(DP-AD)
                oth_locus.append(OTH)

            ref = np.array(ref_locus)
            alt = np.array(alt_locus)
            oth = np.array(oth_locus)
            if allele in SNP_freq:
                freq = SNP_freq[allele]
            else:
                freq = 0 
            
            if keep_locus(ref,alt,oth,minPercentageCellsExpressed=args.minPer):
                REFs.append(ref)
                ALTs.append(alt)
                loci.append(locus)
                if freq>=0.0005:
                    SNPs.append(locus)


df_ref = pd.DataFrame(REFs,index=loci,columns=barcodes)
df_alt =pd.DataFrame(ALTs,index=loci,columns=barcodes)

df_ref.to_csv(os.path.join(args.o,"ref.csv"))
df_alt.to_csv(os.path.join(args.o,"alt.csv"))
if args.freq is not None:
    np.savetxt(os.path.join(args.o,"SNPs.csv"),np.array(SNPs),fmt='%s')
