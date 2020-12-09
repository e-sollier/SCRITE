## Preprocessing

The data preprocessing for this project is different from conventional scRNAseq work-flows, because here we are interested in variants and not in gene expression. For this small project, I did not write a proper pipeline using a tool like nextflow or snakemake, but I instead used simple bash scripts to run the different preprocessing steps. I describe here what I have done, although some parts might have to be adapted depending on which type of data is used.

There are two  main  workflows,  depending  on  whether  we  use  smart-like  data  or  10x  data. For smart-seq data, each cell has a corresponding FASTQ file (or 2 for paired-end sequencing), while for 10x sequencing all of the reads for all the cells are grouped into the same 2 FASTQ files:  the first FASTQ contains the cell barcode and the UMI, while the other FASTQ contains the corresponding reads.

### Tools required
Here are the tools that I used for the preprocessing. Note that this preprocessing is supposed to be done on a cluster, because a large amount of storage is required, and some tools (STAR in particular) require more than 30GB of RAM. Also, each of these tools has to be added to the path after being installed (except cellSNP) by appending `export PATH="<path_to_tool>:$PATH"` to .bashrc.
* [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc). This is used for downloading FASTQ files stored in the Sequence Read Archive. There instructions for installing it are available [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit).

* [STAR](https://github.com/alexdobin/STAR): The [GitHub repository](https://github.com/alexdobin/STAR) contains instructions for the installations and the [manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) explains how to run commands. In particular, it is necessary to first create an index before performing any alignment. You first have to download a reference human genome (for example `wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`) and a gtf file (for example `wget ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz`). Then, the index can be created by running `STAR -runMode genomeGenerate --genomeDir /path/to/genomeDir --genomeFastaFiles /path/to/genome/fasta1 --sjdbGTFfile /path/to/annotations.gtf)`. Be sure to request a lot of memory when submitting the job (ie `bsub -R "rusage[mem=32000]" ...`). Once this is done, STAR is ready to be used for alignment.

* Cellranger (only for 10x data). This is used for aligning 10x data (and also inferring cell barcodes). Instructions for installing cellranger are available [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_in).

* [samtools](http://www.htslib.org/). This is used to manipulate the BAM files. The instructions for installing samtools are available [here](http://www.htslib.org/download/).

* [VarScan](http://varscan.sourceforge.net/). This is used to detect variants. The latest jar file of VarScan can be downloaded from [here](https://sourceforge.net/projects/varscan/files/). Each time you run VarScan, you have to reference this jar file. In addition, to be able to use JAVA on the EULER cluster, you should add the line `module load java/1.8.0_91` to your .bashrc file.

* [cellSNP](https://github.com/single-cell-genetics/cellSNP). This is used to get the ref/alt counts for numerous single cells, when they are in a single BAM file.



### Downloading from SRA

For all of the datasets that I used,  the FASTQs were available at the Sequence ReadArchive (SRA). The [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc) has to be installed. Once  this  is  done,  FASTQs  can  be downloaded with the command `fastq-dump –split-files SRA_id`.  The option `--split-files` ensures that, for paired-end sequencing, both ends are separated into 2 files.  Otherwise, both ends would be concatenated in a single FASTQ file, as if there was no gap between them,  which will confuse the aligner.  For 10x data,  the FASTQ files have to be renamed to the format SampleName_S1_L001_ReadType_001.fastq (where ReadType is (I1), R1 or R2), so that cellranger can recognize them. The naming conventions used by cellranger are described [here](https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-). For smart-seq datasets, it is possible to download an accession list using the SRA run selector. Using the accession list, the script `run_downloadSRA.sh` can be used to download all the FASTQ files.


### Alignment

For smart-seq data, we can directly align the FASTQ files with STAR. For 10x data, we have to use cellranger, which uses STAR under the hood but also assigns reads to cells.

For smart-seq datasets, the script `run_star.sh` can be used to run star for all cells in parallel. It basically amounts to running in parallel the command `STAR --runThreadN 8 --genomeDir ../refdata-gex-GRCh38-2020-A/star_genome --sjdbGTFfile ../refdata-gex-GRCh38-2020-A/genes/genes.gtf --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix alignments/${line} --readFilesIn fastq/${line}_1.fastq fastq/${line}_2.fastq`.

For 10x data, the command looks like: `cellranger count --id out_dir --fastqs fastq --sample SRR12603782 --transcriptome ../GRCh38 --chemistry threeprime` where `fastq` is the directory containing the renamed fastq files.


### Obtaining the count matrices

After the alignment, we get one BAM file per cell with smart-seq data, and one BAM file for all of the cells with 10x data, where in each line a tag indicates the corresponding cell barcode.  For smart-seq data, I merge all the BAM files into one, where I keep the cell information in a read group tag.  That way, both the smart-seq and 10x workflows converge, which make the subsequent analysis easier.  This is done by running `samtools merge -r -b bam_list.txt` where `bam_list.txt` contains the list of the BAM files. In practice, I split the BAM file into several files (one per chromosome) to enable some parallelization. This is done with my script `run_merge.sh`.

Then, I run Varscan on this merged BAM file to identify positions where there might be some variants in some cells.  I filter each locus based on the total coverage and number of alternative reads, but at this point, I am not strict with the filtering. I also use my script `filter_variants.py` to filter out variants which do not have any ref reads (homozygous alt), because VarScan does not allow this type of filter. This is done with the script script `run_varscan.sh`.


Then, I use [cellSNP](https://github.com/single-cell-genetics/cellSNP) to get the counts of reference and alternative reads for each cell and each position that was identified by VarScan. This is done with my script `run_cellSNP.sh`.

The output of cellSNP has to be uncompressed, and then it can be processed by my script `read_cellSNP.py`. It can be run with: `python preprocessing/read_cellSNP.py -i cellSNP_output/glioblastoma_BT_S2 -o Data/glioblastoma_BT_S2 --freq genome1K`. genome1K is a directory containing one VCF file per chromosome, with the frequency of the common SNPs in the population. It was created by downloading the file genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz from https://sourceforge.net/projects/cellsnp/files/SNPlist/ and then running my script `split_genome1K.py`. It can be run with: `python split_genome1K -i path/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf -o ../genome1K`.
