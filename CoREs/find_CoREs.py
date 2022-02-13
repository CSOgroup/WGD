#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script to identify Compartment repositioning events from pairs of Calder 
segmentations in sub-compartments.


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
import numpy as np
import pandas as pd
from pybedtools.bedtool import BedTool


CALDER_SUBCOMPARTMENT_HEADER = ['chr', 'start', 'end', 'compartment_label', 'domain_rank_8_levels', 
								'strand', 'tickStart', 'tickEnd', 'color']


def load_rank_full(sample_path, binsize):
	"""Given a path to a Calder generated segmentation file (.bed format), it loads it and
	returns a bin-level table where each bin is associated to its Calder rank.
	It accepts:
		- sample_path: path to the Calder generated .bed file
		- binsize: resolution at which sub-compartments were called from Hi-C
	"""

	s_comps = pd.read_csv(sample_path, sep='\t', header=None, names = CALDER_SUBCOMPARTMENT_HEADER)
	s_comps['compartment_label_8'] = s_comps.compartment_label.map(lambda x: x[:5])

	domain_rank = s_comps.groupby("chr").compartment_label.rank(ascending=False, pct=True)
	s_comps['domain_rank'] = domain_rank

	min_max_domain_rank = s_comps.groupby('chr')['domain_rank']\
								 .min()\
								 .to_frame("min_domain_rank")\
								 .reset_index()\
								 .merge(s_comps.groupby('chr')['domain_rank'].max().to_frame("max_domain_rank").reset_index())

	s_comps = s_comps.merge(min_max_domain_rank)
	s_comps['domain_rank'] = (s_comps.domain_rank - s_comps.min_domain_rank)/(s_comps.max_domain_rank - s_comps.min_domain_rank)

	s_rank = s_comps[['chr','start','end', 'domain_rank']]
	bins = BedTool().window_maker(w=binsize, genome='hg19')
	s_rank = bins.map(BedTool.from_dataframe(s_rank).sort(), c=4, o='max', null=np.NaN)\
				 .to_dataframe(names = ['chr','start','end', 'domain_rank'])
	s_rank = s_rank[s_rank.chr.isin(s_comps.chr.unique())].reset_index(drop=True)
	return s_rank




def recursive_segmentation(df, 
						   min_std = 0.1, 
						   control_dist=None, 
						   value_col = 'delta_rank', 
						   **kwrgs):

	"""Implements the recursive segmentation strategy which, given a genomic signal representing the delta compartment rank
	for each genomic bin, splits the genome in regions on the basis of their average delta-rank.
	Each segment is associated to its average delta rank. Segments are determined on the basis of their standard deviation, which
	cannot exceed a pre-defined threshold. If a control distribution of delta ranks is provided, an empirical p-value is also calculated
	for each detected segment. This function accepts:
		- df: pandas dataframe with columns chr, start, end, delta-rank. It must be ordered and bins have to have the same size and being adjacent
		- min_std: minimum standard deviation required for splitting a segment
		- control_dist [optional]: dictionary having chromosomes as keys and vectors of delta-ranks as values, representing
									the expected values calculated from control comparisons
		- value_col: column of df corresponding to delta-rank
		- kwargs: other parameters passed to the function _pos_neg_segmentation
	"""

	control_dist_abs = {chrom:np.abs(v) for chrom, v in control_dist.items()} if control_dist is not None else None


	def __pos_neg_segmentation(v, segment_value_function = np.mean):
		""" Given a signal, it segments it based on positive and negative values
		It accepts:
		- v: Numpy array
		- segment_value_function: function accepting a numpy array and returning the imputed
			value for each segment, on the basis of the array values
		- 
		"""

		segments = []
		seg_values = []
		seg_start = 0
		prev_v = np.NaN
		for i in range(v.shape[0]):
			if not np.isnan(prev_v):
				if np.sign(v[i]) != np.sign(prev_v):
					segments.append([seg_start, i])
					seg_values.append(segment_value_function(v[seg_start:i]))
					seg_start = i
			else:
				seg_start = i
			prev_v = v[i]
		if not np.isnan(prev_v):
			segments.append([seg_start, i + 1])
			seg_values.append(segment_value_function(v[seg_start:i + 1]))
		segments = np.array(segments, dtype=int)
		seg_values = np.array(seg_values)
		return segments, seg_values


	def __recursive_segmentation(X, sub_mean=True):
		X = X.sort_values(['chr', 'start', 'end']).reset_index(drop=True)
		v = X[value_col].values
		mean_signal = np.nanmean(v)
		std_signal = np.nanstd(v)
		if std_signal > min_std:
			segments, _ = __pos_neg_segmentation(v - mean_signal if sub_mean else v, **kwrgs)
			all_rrs = []
			for s in range(segments.shape[0]):
				seg = segments[s, :]
				Xr = X.iloc[seg[0]:seg[1]]
				rrs = __recursive_segmentation(Xr)
				all_rrs.append(rrs)
			result = pd.concat(all_rrs, axis=0, ignore_index=True)
		else:
			result = BedTool.from_dataframe(X.dropna(subset=[value_col])).merge().to_dataframe(names = ['chr', 'start', 'end'])
			result['value'] = mean_signal

			if (control_dist is not None) and ((~np.isnan(v)).sum() > 0):
				pvalue = ((control_dist_abs[X.chr.iloc[0]] > np.max(np.abs( v[~np.isnan(v)] ))).sum() + 1)/(control_dist_abs[X.chr.iloc[0]].shape[0] + 1)
				result['empirical_pvalue'] = pvalue
			else:
				result['empirical_pvalue'] = np.NaN
		result = result.sort_values(['chr', 'start', 'end'])
		return result


	result = []
	for chrom, X in df.groupby('chr'):
		chrom_result = __recursive_segmentation(X, sub_mean=False)
		result.append(chrom_result)
	result = pd.concat(result, axis=0, ignore_index=True)
	if control_dist is None:
		result = result.drop("empirical_pvalue", axis=1)
	return result


def create_parser():
	parser = argparse.ArgumentParser(description = "Identifying Compartment Repositioning Events from Calder genomic segmentations")
	parser.add_argument("sample1_path", type=str, help="Path to the Calder segmentation of sample 1")
	parser.add_argument("sample2_path", type=str, help="Path to the Calder segmentation of sample 2")
	parser.add_argument("binsize", type=int, help="Resolution to use in the analysis")
	parser.add_argument("output_path", type=str, help="Path where to store the identified regions")
	parser.add_argument("--min_std", type=float, default=0.1, help="Maximum standard deviation allowed for segmented regions")
	parser.add_argument("--control1_path", type=str, nargs='*', help="Path(s) to the Calder segmentation(s) to use to use as control 1")
	parser.add_argument("--control2_path", type=str, nargs='*', help="Path(s) to the Calder segmentation(s) to use to use as control 2")
	return parser


def main():
	parser = create_parser()
	args = parser.parse_args()

	null_dist = None
	if (args.control1_path is not None) and (args.control2_path is not None):
		all_controls = []
		for c1_path, c2_path in zip(args.control1_path, args.control2_path):
			c1_rank = load_rank_full(c1_path, args.binsize)
			c2_rank = load_rank_full(c2_path, args.binsize)
			c12_rank = c1_rank.merge(c2_rank, on=['chr', 'start', 'end'], suffixes=("_1", "_2"))
			c12_rank['delta_rank'] = c12_rank.domain_rank_2 - c12_rank.domain_rank_1
			c12_rank['comparison'] = f"{c1_path} VS {c2_path}"
			all_controls.append(c12_rank)
		all_controls = pd.concat(all_controls, axis=0, ignore_index=True)
		null_dist = all_controls.dropna(subset=['delta_rank']).groupby('chr').apply(lambda x: x.delta_rank.values).to_dict()


	s1_rank = load_rank_full(args.sample1_path, args.binsize)
	s2_rank = load_rank_full(args.sample2_path, args.binsize)
	s12_rank = s1_rank.merge(s2_rank, on=['chr', 'start', 'end'], suffixes=("_1", "_2"))
	s12_rank['delta_rank'] = s12_rank.domain_rank_2 - s12_rank.domain_rank_1
	segmentation = recursive_segmentation(s12_rank, min_std=args.min_std, value_col='delta_rank', control_dist = null_dist)
	segmentation.to_csv(args.output_path, sep="\t", index=False, header=True)



if __name__ == '__main__':
	main()







