# SCRITE

SCRITE (Single-Cell Rna-seq Inference of Tumour Evolution) is an algorithm which estimates the phylogeny of the somatic events (mutations and LOH) in a tumor, using scRNAseq data. 
It is based on [SCITE](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0936-x) / [SCIPhi](https://www.nature.com/articles/s41467-018-07627-7). The Python implementation was adapted from https://github.com/znorio/RNAmut.


## Dependencies

The file SCRITE.yml contains the dependencies. They can be installed by running `conda env create -f SCRITE.yml`.

## Demo

The jupyter notebook demo.ipynb shows an example on one dataset.

## Input Data

SCRITE takes as input:
* ref.csv: matrix of the counts for the reference allele (rows: mutations, columns: cells)
* alt.csv: matrix of the counts for the alternate allele (rows: mutations, columns: cells)
* SNPs.csv (optional): list of the mutations present in ref.csv and alt.csv which are SNPs (i.e. have a high frequency in the population). This makes the selection of candidate LOH events and somatic mutations easier.
* metadata.csv (optional): metadata for the cells, for example the true label (Neoplastic/Regular) to compare the results

In this repository, the Data directory contains one example preprocessed dataset. The other preprocessed datasets that I used are available at: https://polybox.ethz.ch/index.php/s/bXkQIh25rq9lqgQ.

The preprocessing directory contains instructions to preprocess a dataset, starting from FASTQ files.
