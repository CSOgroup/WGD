#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Aggregating Hi-C interactions by bin categories.
   This script takes as input:
   - an Hi-C matrix path, in Cooler format
   - the path to a bin file, in tsv format, having as 
   	 columns [chr, start, end, category], without header,
   	 where each line is a Hi-C bin together with their category
   - the output path where to store the result, which is a 
     tsv file, having as columns [chr1, category1, chr2, category2, n_interactions]
"""



__author__ 		=	"Luca Nanni"
__contact__ 	=	"luca.nanni@unil.ch"
__date__ 		=	"2022/02/13"
__version__ 	=	"0.0.1"



import pandas as pd
from pybedtools.bedtool import BedTool
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
from cooler import Cooler
import argparse


CALDER_SUBCOMPARTMENT_HEADER = ['chr', 'start', 'end', 'compartment_label', 'domain_rank_8_levels', 
								'strand', 'tickStart', 'tickEnd', 'color']


def process_arguments():
	parser = argparse.ArgumentParser(description='Coarse Hi-C data by bin categories')
	parser.add_argument("cool_path", type=str, help='Path to the cool file')
	parser.add_argument("resolution", type=int, help='Resolution to use for the Hi-C file')
	parser.add_argument("compartments_path", type=str, help='Path to the bin category file')
	parser.add_argument('output_path', type=str, help='Output path where to store the category aggregated result')
	args = parser.parse_args()
	return args


def load_compartments(compartments_path, binsize):
	s_comps = pd.read_csv(compartments_path, sep='\t', header=None, names = CALDER_SUBCOMPARTMENT_HEADER)
	chroms = s_comps.chr.unique()
	s_comps['category'] = s_comps.compartment_label.map(lambda x: x[:5])
	s_comps = s_comps[['chr', 'start', 'end', 'category']]
	bins = BedTool().window_maker(w=binsize, genome='hg19')
	s_comps = bins.map(BedTool.from_dataframe(s_comps).sort(), c=4, o='distinct', null="Unknown")\
				 .to_dataframe(names = ['chr','start','end', 'category'])
	s_comps = s_comps[s_comps.chr.isin(chroms)].reset_index(drop=True)
	return s_comps


def aggregate(cool_ints, bins, categories):
	def __aggregate_intrachromosomal(cool_ints):
		cool_cat_aggregation = cool_ints.groupby(['chr1', 'chr2', 'category1', 'category2'])['balanced'].sum().compute()
		cool_cat_aggregation = cool_cat_aggregation.to_frame('count').reset_index()
		cool_cat_aggregation['category1'] = pd.Categorical(cool_cat_aggregation['category1'], categories=categories, ordered=True)
		cool_cat_aggregation['category2'] = pd.Categorical(cool_cat_aggregation['category2'], categories=categories, ordered=True)
		cat1 = cool_cat_aggregation.apply(lambda x: min(x.category1, x.category2), axis=1)
		cat2 = cool_cat_aggregation.apply(lambda x: max(x.category1, x.category2), axis=1)
		cool_cat_aggregation['category1'] = cat1
		cool_cat_aggregation['category2'] = cat2
		cool_cat_aggregation = cool_cat_aggregation.groupby(['chr1', 'chr2', 'category1', 'category2'])['count'].sum()
		cool_cat_aggregation = cool_cat_aggregation.to_frame('count').reset_index()
		return cool_cat_aggregation

	def __aggregate_interchromosomal(cool_ints):
		cool_cat_aggregation = cool_ints.groupby(['chr1', 'chr2', 'category1', 'category2'])['balanced'].sum().compute()
		cool_cat_aggregation = cool_cat_aggregation.to_frame('count').reset_index()
		return cool_cat_aggregation


	cool_ints = cool_ints.merge(bins[['chr','start','category']].rename(columns=lambda x: x + "1"), on=['chr1', 'start1'])
	cool_ints = cool_ints.merge(bins[['chr','start','category']].rename(columns=lambda x: x + "2"), on=['chr2', 'start2'])

	cool_cat_aggregation = pd.concat([
		__aggregate_intrachromosomal(cool_ints[cool_ints.chr1 == cool_ints.chr2]), 
		__aggregate_interchromosomal(cool_ints[cool_ints.chr1 != cool_ints.chr2]), 
		], axis=0, ignore_index=True)

	cool_cat_aggregation = cool_cat_aggregation[['chr1', 'category1', 'chr2', 'category2', 'count']]
	return cool_cat_aggregation


def aggregate_HiC_counts_by_category(cool, bins):
	bins = dd.from_pandas(bins, npartitions=10)
	categories = sorted(bins.category.unique().compute())

	cool_ints = cool.matrix(as_pixels=True, join=True, balance='KR')[:]
	cool_ints['chrom1'] = cool_ints['chrom1'].astype(str)
	cool_ints['chrom2'] = cool_ints['chrom2'].astype(str)

	x = cool_ints.iloc[0].chrom1
	
	cool_ints = dd.from_pandas(cool_ints, npartitions=10)
	
	if not x.startswith('chr'):
		cool_ints['chrom1'] = 'chr' + cool_ints.chrom1
		cool_ints['chrom2'] = 'chr' + cool_ints.chrom2
	
	cool_ints = cool_ints.rename(columns={'chrom1':'chr1', 'chrom2': 'chr2'})
	cool_cat_aggregation = aggregate(cool_ints, bins, categories)
	return cool_cat_aggregation

def main():
	ProgressBar().register()
	args = process_arguments()

	compartments = load_compartments(args.compartments_path, args.resolution)
	cool = Cooler(args.cool_path + f"::/resolutions/{args.resolution}")
	cool_cat_aggregation = aggregate_HiC_counts_by_category(cool, compartments)
	cool_cat_aggregation.to_csv(args.output_path, sep='\t', header=True, index=False)

if __name__ == '__main__':
	main()














