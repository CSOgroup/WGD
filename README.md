# WGD

Analysis code for the paper 

> *Whole genome doubling promotes tumorigenesis through loss of chromatin segregation and compartment repositioning* 
> 
> by Ruxandra A. Lambuta, Luca Nanni, Yuanlong Liu, Juan Diaz-Miyar, Arvind Iyer, Daniele Tavernari, Natalya Katanayeva, Giovanni Ciriello and Elisa Oricchio


## Hi-C processing

### Converting `.hic` to `.cool`
Most of the Hi-C analyses done in the paper use the `.cool` format. To convert the `.hic` files stored on Zenodo ([link1](https://zenodo.org/record/6053792), [link2](https://zenodo.org/record/6054423)), you can use the following script:

```
bash hic_processing/convert_hic_to_cool.sh input.hic output.mcool
```

## Interactions between chromosomes at WGD
Given a Hi-C experiment in `.cool` format, we can dump interactions at the chromosome level using the following script:
```
python interchromosomal/dump_cool.py --resolution 1000000000 --genome hg19 input.cool interchromosomal_dump.txt
```
The script will store the inter-chromosomal interactions in a text file and also perform IC balancing of the contacts.

Given two interchromosomal contact files `sample1_dump.txt` and `sample2_dump.txt`, we can now compare them using:
```
python interchromosomal/compare_interchromosomal_contacts.py \
                               sample1_dump.txt \
                               sample2_dump.txt \
                               interchromosomal_output \
                               --name1 sample1 \
                               --name2 sample2
```

Explanation of the parameters:
```
usage: compare_interchromosomal_contacts.py [-h] [--name1 NAME1] [--name2 NAME2] sample1 sample2 output_path

Compares the inter-chromosomal contacts between two Hi-C experiments

positional arguments:
  sample1        Path to the first sample (control)
  sample2        Path to the second sample (treatment)
  output_path    Where to store the results

optional arguments:
  -h, --help     show this help message and exit
  --name1 NAME1  Name of first sample
  --name2 NAME2  Name of second sample
```

This command will create a folder with four plots inside:
* Observed/Expected inter-chromosomal contacts for sample1
* Observed/Expected inter-chromosomal contacts for sample2
* Contact fold-change between sample1 and sample2 (sample2  / sample1)
* Boxplot comparing the fold-change distributions between long and short chromosome pairs



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