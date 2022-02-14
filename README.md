# WGD

Analysis code for the paper 

> *Whole genome doubling promotes tumorigenesis through loss of chromatin segregation and compartment repositioning* 
> 
> by Ruxandra A. Lambuta, Luca Nanni, Yuanlong Liu, Juan Diaz-Miyar, Arvind Iyer, Daniele Tavernari, Natalya Katanayeva, Giovanni Ciriello and Elisa Oricchio

## Calder compartment calls
All the Hi-C sub-compartment calls in this study made with Calder are deposited at [this Zenodo link](https://zenodo.org/record/6054423). BED files used in the analysis are in the subfolder `hic_features/compartment_domains`.

## Detecting Compartment Repositioning Events (CoREs)
Compartment repositioning events are detected directly from a pair of Calder compartment call files (.bed).

```
python CoREs/find_CoREs.py sample1_calder.bed \
			   sample2_calder.bed \
			   50000 \
			   sample1_sample2_CoREs.tsv \
			   --min_std 0.1 \
			   --control1_path control1_calder.bed \
			   --control2_path control2_calder.bed
```

Meaning of the arguments:

```
find_CoREs.py [-h] [--min_std MIN_STD] [--control1_path [CONTROL1_PATH ...]] [--control2_path [CONTROL2_PATH ...]] sample1_path sample2_path binsize output_path

Identifying Compartment Repositioning Events from Calder genomic segmentations

positional arguments:
  sample1_path          Path to the Calder segmentation of sample 1
  sample2_path          Path to the Calder segmentation of sample 2
  binsize               Resolution to use in the analysis
  output_path           Path where to store the identified regions

optional arguments:
  -h, --help            show this help message and exit
  --min_std MIN_STD     Maximum standard deviation allowed for segmented regions
  --control1_path [CONTROL1_PATH ...]
                        Path(s) to the Calder segmentation(s) to use to use as control 1
  --control2_path [CONTROL2_PATH ...]
                        Path(s) to the Calder segmentation(s) to use to use as control 2
```