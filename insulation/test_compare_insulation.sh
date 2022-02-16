#!/bin/bash


sample1_insulation="../wgd/results/hic/insulation/RPE_Control_insulation_w1000000_r50000_KR.bw"
sample1_boundaries='../wgd/results/hic/insulation/RPE_Control_insulation_w1000000_r50000_KR_boundaries.bed'
sample2_insulation="../wgd/results/hic/insulation/RPE_Treated1_insulation_w1000000_r50000_KR.bw"
sample2_boundaries='../wgd/results/hic/insulation/RPE_Treated1_insulation_w1000000_r50000_KR_boundaries.bed'
output_path="./out"
name1="RPE_Control"
name2="RPE_WGD1"

python insulation/compare_insulation.py ${sample1_insulation} \
										${sample1_boundaries} \
										${sample2_insulation} \
										${sample2_boundaries} \
										${output_path} \
										--name1 ${name1} \
										--name2 ${name2}