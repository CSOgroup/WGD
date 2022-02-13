#!/bin/bash

sample1_path="../wgd/data/published/zenodo/hic/hic_features/compartment_domains/RPE_TP53_Ctrl_mega_CompartmentDomains.bed"
sample2_path="../wgd/data/published/zenodo/hic/hic_features/compartment_domains/RPE_TP53_20w0T1_mega_CompartmentDomains.bed"
binsize=50000
output_path="./out.tsv"

control1_path="../wgd/data/published/zenodo/hic/hic_features/compartment_domains/RPE_TP53_Ctrl_1_CompartmentDomains.bed"
control2_path="../wgd/data/published/zenodo/hic/hic_features/compartment_domains/RPE_TP53_Ctrl_2_CompartmentDomains.bed"
min_std=0.1


python CoREs/find_CoREs.py ${sample1_path} \
						   ${sample2_path} \
						   ${binsize} ${output_path} \
						   --min_std ${min_std} \
						   --control1_path ${control1_path} \
						   --control2_path ${control2_path}