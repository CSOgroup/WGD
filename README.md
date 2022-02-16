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


## Interactions between sub-compartments at WGD

Given a Hi-C experiment `input.mcool` and its relative Calder sub-compartments `input_calder.bed`, we can then aggregate the contacts for each pair of Calder sub-compartments levels, for each chromosome and chromosome pair as follows:

```
python intercompartments/aggregate_interactions_by_compartment.py \
        input.mcool \
        50000 \
        input_calder.bed \
        input_interactions_by_compartment.tsv
```

Meaning of the arguments:
```
usage: aggregate_interactions_by_compartment.py [-h] cool_path resolution compartments_path output_path

Coarse Hi-C data by bin categories

positional arguments:
  cool_path          Path to the cool file
  resolution         Resolution to use for the Hi-C file
  compartments_path  Path to the bin category file
  output_path        Output path where to store the category aggregated result

optional arguments:
  -h, --help         show this help message and exit
```

Once aggregated the interactions at the compartment level, we can compare two samples `sample1_interactions_by_compartment.tsv` and `sample2_interactions_by_compartment.tsv` as follows:

```
python intercompartments/compare_intercompartments_contacts.py \
        sample1_interactions_by_compartment.tsv \
        sample2_interactions_by_compartment.tsv \
        intercompartments_output \
        --name1 sample1 \
        --name2 sample2
```

Meaning of the arguments:
```
usage: compare_intercompartments_contacts.py [-h] [--name1 NAME1] [--name2 NAME2] sample1 sample2 output_path

Coarse Hi-C data by bin categories

positional arguments:
  sample1        Path to compartments contacts for the first sample (control)
  sample2        Path to compartments contacts for the first sample (treatment)
  output_path    Where to store the results

optional arguments:
  -h, --help     show this help message and exit
  --name1 NAME1  Name of first sample
  --name2 NAME2  Name of second sample
```

This command will create a folder with six plots inside:
* Observed/Expected inter-compartments contacts for sample1 focusing on intra-chromosomal interactions
* Observed/Expected inter-compartments contacts for sample1 focusing on inter-chromosomal interactions
* Observed/Expected inter-compartments contacts for sample2 focusing on intra-chromosomal interactions
* Observed/Expected inter-compartments contacts for sample2 focusing on inter-chromosomal interactions
* Inter-compartments contact fold-change between sample1 and sample2 (sample2  / sample1) focusing on intra-chromosomal interactions
* Inter-compartments contact fold-change between sample1 and sample2 (sample2  / sample1) focusing on inter-chromosomal interactions


## Hi-C Insulation scores and insulation boundaries
All the Hi-C insulation scores for each sample of this study are deposited at [this Zenodo link](https://zenodo.org/record/6054423). BED files used in the analysis are in the subfolder `hic_features/insulation_scores`.

Hi-C boundaries derived from insulation scores are available in the subfolder `hic_features/insulation_boundaries`.

## Comparing Hi-C insulation between Control and WGD
Given the insulation scores and the boundairies of two samples, we can compare their levels across two conditions as follows:
```
python insulation/compare_insulation.py \
                    sample1_insulation.bw \
                    sample1_boundaries.bed \
                    sample2_insulation.bw \
                    sample2_boundaries.bed \
                    insulation_output \
                    --name1 sample1 \
                    --name2 sample2
```

Meaning of the arguments:
```
usage: compare_insulation.py [-h] [--name1 NAME1] [--name2 NAME2]
                             sample1_insulation sample1_boundaries sample2_insulation sample2_boundaries output_path

Compares the insulation between two samples

positional arguments:
  sample1_insulation  Path to the first sample insulation (control)
  sample1_boundaries  Path to the first sample boundaries (control)
  sample2_insulation  Path to the second sample (treatment)
  sample2_boundaries  Path to the second sample boundaries (treatment)
  output_path         Where to store the results

optional arguments:
  -h, --help          show this help message and exit
  --name1 NAME1       Name of first sample
  --name2 NAME2       Name of second sample
```

This command will output one plot:
* Scatterplot where for each shared boundary between the two conditions the insulation values are shown


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