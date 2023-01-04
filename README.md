# Whole genome doubling drives oncogenic loss of chromatin segregation

Ruxandra A. Lambuta, Luca Nanni, Yuanlong Liu, Juan Diaz-Miyar, Arvind Iyer, Daniele Tavernari, Natalya Katanayeva, Giovanni Ciriello and Elisa Oricchio


# Running the analysis
You can run the whole set of analyses by simply running 
```
bash run.sh
```

The script will:
1. Download the data necessary for the analysis from Zenodo
2. Run each analysis present in the [`src`](./src) folder 
3. Store the resulting figures in the [`figures`](./figures) folder

You can check the rest of this document to learn about each step of the pipeline.

## [00_download.sh](src/00_download.sh)
Downloads the data necessary for the analysis from the two Zenodo repositories:
- Part 1: https://doi.org/10.5281/zenodo.7351767
- Part 2: https://doi.org/10.5281/zenodo.7351776

It then merges the files in a single `zip` file, and finally decompress the archive into the `data` folder

## [01_hic_to_mcool.sh](src/01_hic_to_mcool.sh)
Converts the `*.hic` files storing information about Hi-C experiments into `*.mcool` files for later processing.

## [02_LCS.py](src/02_LCS.py)
It performs, for each sample comparison presented in the paper, the analyis of Loss of Chromatin Segregation (LCS).
The analysis is done a three different levels:
- Inter-chromosomal interactions
- Inter-compartmental interactions
- Boundary insulation

## [03_ChIPSeq.py](src/03_ChIPSeq.py)
It performs the analysis of ChIPSeq signals and peaks for CTCF and H3K9me3 between Control and WGD samples, as presented in the paper. 

## [04_scHiC.py](04_scHiC.py)
Analysis of Single-cell Hi-C data, as presented in the paper. Specifically, the code will:
1. Analyse inter-chromosomal interactions in single cells
2. Perform pseudo-bulk inter-chromosomal analysis from single cells
3. Analyse compartment segregation at single cell level
4. Infer copy number changes for each single cell