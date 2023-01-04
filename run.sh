#!/bin/bash


# Downloading Zenodo repository
bash src/00_download.sh

# Converting Hi-C files to MCOOL files for later usage
bash src/01_hic_to_mcool.sh

# Computing LCS measurements on the Hi-C samples
python src/02_LCS.pys

# Analysing ChIP-seq data
python src/03_ChIPSeq.py

# Analysing Sc-HiC data
python src/04_scHiC.py