#!/bin/bash


sample1_path="../wgd/results/hic_compartments_aggregations/RPE_Control_mindist_0.txt"
sample2_path="../wgd/results/hic_compartments_aggregations/RPE_Treated1_mindist_0.txt"
output_path="./out"
name1="RPE_Control"
name2='RPE_WGD1'


python intercompartments/compare_intercompartments_contacts.py \
				${sample1_path} \
				${sample2_path} \
				${output_path} \
				--name1 ${name1} \
				--name2 ${name2}