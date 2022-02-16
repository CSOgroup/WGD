#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Compares Hi-C insulation and insulation boundaries between two samples

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details
"""

__author__ 		=	"Luca Nanni"
__contact__ 	=	"luca.nanni@unil.ch"
__date__ 		=	"2022/02/13"
__version__ 	=	"0.0.1"


import argparse
import os
import pyBigWig
import pandas as pd
from pybedtools.bedtool import BedTool
import matplotlib.pyplot as plt
import seaborn as sns



CHROMS = ['chr' + str(x) for x in range(1, 23)] + ['chrX']
TOP_K=300
insulation_change_to_color = {
		'gain': 'red',
		'loss': 'blue'
	}


def process_arguments():
    parser = argparse.ArgumentParser(description='Compares the insulation between two samples')
    parser.add_argument("sample1_insulation", type=str, help='Path to the first sample insulation (control)')
    parser.add_argument("sample1_boundaries", type=str, help='Path to the first sample boundaries (control)')
    parser.add_argument('sample2_insulation', type=str, help='Path to the second sample (treatment)')
    parser.add_argument("sample2_boundaries", type=str, help='Path to the second sample boundaries (treatment)')
    parser.add_argument("output_path", type=str, help="Where to store the results")
    parser.add_argument("--name1", type=str, default="sample1", help="Name of first sample")
    parser.add_argument("--name2", type=str, default="sample2", help="Name of second sample")
    args = parser.parse_args()
    return args


def load_bigwig(path):
	bw = pyBigWig.open(path)

	intervals = []
	for chrom in bw.chroms():
		chrom_intervals = bw.intervals(chrom)
		chrom_intervals = pd.DataFrame.from_records(chrom_intervals, columns = ["start", "end", "score"])
		if not chrom.startswith("chr"):
			chrom = "chr" + chrom
		chrom_intervals['chr'] = chrom
		chrom_intervals = chrom_intervals[['chr', 'start', 'end', 'score']]
		intervals.append(chrom_intervals)
	intervals = pd.concat(intervals, axis=0, ignore_index=True)
	intervals = intervals[intervals.chr.isin(CHROMS)].reset_index(drop=True)
	return intervals


def load_boundaries(path):
	sample_boundaries = pd.read_csv(path, sep="\t", header=None, names = ['chr', 'start', 'end', 'name', 'score', 'strand'])
	sample_boundaries['chr'] = sample_boundaries.chr.astype(str)
	if not sample_boundaries.chr.iloc[0].startswith("chr"):
		sample_boundaries['chr'] = 'chr' + sample_boundaries.chr.astype(str)
	sample_boundaries = sample_boundaries[sample_boundaries.chr.isin(CHROMS)]
	sample_boundaries['start'] -= 1
	sample_boundaries = sample_boundaries[['chr', 'start', 'end', 'score']]
	return sample_boundaries


def get_shared_boundaries(control_boundaries, control_insulation, wgd_boundaries, wgd_insulation):
	shared = BedTool.from_dataframe(control_boundaries[['chr', 'start', 'end', 'score']].sort_values(['chr', 'start', 'end']))\
		   			.closest(BedTool.from_dataframe(wgd_boundaries[['chr', 'start', 'end', 'score']].sort_values(['chr', 'start', 'end'])), d=True)\
		   			.to_dataframe(names = ['chr', 'start', 'end', 'control_score', "wgd_chr", "wgd_start", "wgd_end", "wgd_score", "distance"])
	shared = shared[shared.distance <= 1]
	shared = shared.merge(control_insulation.rename(columns = {'score': 'control_ins'}), on=['chr', 'start', 'end'])\
		  		   .merge(wgd_insulation.rename(columns = {'score': 'wgd_ins'}), on=['chr', 'start', 'end'])
	shared['insulation_change'] = shared.apply(lambda y: "loss" if y.wgd_ins > y.control_ins else "gain", axis=1)
	shared['control_ins_rank'] = shared.control_ins.rank()
	shared['wgd_ins_rank'] = shared.wgd_ins.rank()
	return shared


def plot_shared_boundaries_insulation(shared_boundaries, output_path):
	with sns.plotting_context('paper', font_scale = 1.3):
		fig = plt.figure()
		sns.scatterplot(data=shared_boundaries, 
					x='control_ins', 
					y='wgd_ins', 
					hue='insulation_change', 
					palette = insulation_change_to_color,
					alpha=0.4,
					rasterized=False)
		xmin, xmax = min(shared_boundaries.control_ins.min(), shared_boundaries.wgd_ins.min()), max(shared_boundaries.control_ins.max(), shared_boundaries.wgd_ins.max())
		plt.plot([xmin, xmax], [xmin, xmax], linestyle = '--', color = 'black', linewidth=2)
		df = shared_boundaries[(shared_boundaries.control_ins_rank <= TOP_K) & (shared_boundaries.wgd_ins_rank <= TOP_K)]
		# find top-right corner
		v = max(df.control_ins.max(), df.wgd_ins.max())
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()
		plt.vlines([v], ymin, v, color = 'black', linestyle = '-', linewidth=1)
		plt.hlines([v], xmin, v, color = 'black', linestyle = '-', linewidth=1)
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.text(x=v - 0.8, y = ymin + 0.4, s = f"Top-{TOP_K}\nboundaries", ha='center')

		plt.xlabel("Control insulation")
		plt.ylabel("WGD insulation")
		sns.despine()
		fig.savefig(output_path, bbox_inches='tight', transparent=True)
		plt.close(fig)



def main():
	args = process_arguments()

	os.makedirs(args.output_path, exist_ok=True)

	sample1_insulation = load_bigwig(args.sample1_insulation).assign(sample = args.name1)
	sample1_boundaries = load_boundaries(args.sample1_boundaries).assign(sample = args.name1)
	sample2_insulation = load_bigwig(args.sample2_insulation).assign(sample = args.name2)
	sample2_boundaries = load_boundaries(args.sample2_boundaries).assign(sample = args.name2)

	shared = get_shared_boundaries(sample1_boundaries, sample1_insulation, sample2_boundaries, sample2_insulation)
	plot_shared_boundaries_insulation(shared, os.path.join(args.output_path, f"{args.name1}_VS_{args.name2}_boundary_insulation.pdf"))




if __name__ == '__main__':
	main()







