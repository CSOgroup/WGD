#!/bin/bash


cool_path="../wgd/results/cool_maps/RPE_Control.mcool"
resolution=50000
compartments_path="../wgd/data/published/zenodo/hic/hic_features/compartment_domains/RPE_TP53_Ctrl_mega_CompartmentDomains.bed"
output_path="./RPE_Control_compartments_interactions.tsv"

python intercompartments/aggregate_interactions_by_compartment.py \
			${cool_path} \
			${resolution} \
			${compartments_path} \
			${output_path}