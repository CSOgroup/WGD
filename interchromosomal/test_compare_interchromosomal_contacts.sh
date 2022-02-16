#!/bin/bash


sample1_path="../wgd/results/hic_dumped_matrices/RPE_Control/dump_1000000000.tsv"
sample2_path="../wgd/results/hic_dumped_matrices/RPE_Treated1/dump_1000000000.tsv"
output_path="./out"
name1="RPE_Control"
name2="RPE_WGD1"

python interchromosomal/compare_interchromosomal_contacts.py ${sample1_path} \
															 ${sample2_path} \
															 ${output_path} \
															 --name1 ${name1} \
															 --name2 ${name2}