#!/usr/bin/env python
# -*- coding: utf-8 -*-


__author__ 		=	"Luca Nanni"
__contact__ 	=	"luca.nanni@unil.ch"
__date__ 		=	"2022/02/13"
__version__ 	=	"0.0.1"


import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


COMPARTMENTS_8_LABELS = ["{}.{}.{}".format(l1, l2, l3) for l1 in ['A', 'B'] for l2 in [1, 2] for l3 in [1, 2]] + ['Unknown']
CHROMS = ['chr{}'.format(x) for x in range(1, 23)] + ['chrX']


def process_arguments():
	parser = argparse.ArgumentParser(description='Coarse Hi-C data by bin categories')
	parser.add_argument("sample1", type=str, help='Path to compartments contacts for the first sample (control)')
	parser.add_argument("sample2", type=str, help='Path to compartments contacts for the first sample (treatment)')
	parser.add_argument("output_path", type=str, help="Where to store the results")
	parser.add_argument("--name1", type=str, default="sample1", help="Name of first sample")
	parser.add_argument("--name2", type=str, default="sample2", help="Name of second sample")
	args = parser.parse_args()
	return args



def intrachromosomal_oe(sample_compartment_agg):
	c_sample_compartment_agg = sample_compartment_agg[sample_compartment_agg.chr1 == sample_compartment_agg.chr2].copy()
	comp1 = c_sample_compartment_agg.compartment1.copy()
	comp2 = c_sample_compartment_agg.compartment2.copy()
	c_sample_compartment_agg_other = c_sample_compartment_agg[c_sample_compartment_agg.compartment1 != c_sample_compartment_agg.compartment2].copy()
	c_sample_compartment_agg_other['compartment1'] = comp2
	c_sample_compartment_agg_other['compartment2'] = comp1
	c_sample_compartment_agg_both = pd.concat([c_sample_compartment_agg, c_sample_compartment_agg_other], axis=0, ignore_index=True)

	c_sample_comp_coverage = c_sample_compartment_agg_both.groupby(['chr1','compartment1', 'chr2'], observed=True)['count']\
								.sum().to_frame('coverage').reset_index().rename(columns={'compartment1': 'compartment'})\
								[['compartment', 'chr1', 'chr2', 'coverage']]
	c_sample_tot_ints = c_sample_compartment_agg.groupby(['chr1', 'chr2'], observed=True)['count'].sum().to_frame('total_coverage').reset_index()
	c_sample_compartment_agg_oe = c_sample_compartment_agg\
				.merge(c_sample_comp_coverage.rename(columns=lambda x: x + "1" if x not in c_sample_compartment_agg.columns else x), on=['compartment1', 'chr1', 'chr2'])\
				.merge(c_sample_comp_coverage.rename(columns=lambda x: x + "2" if x not in c_sample_compartment_agg.columns else x), on=['compartment2', 'chr1', 'chr2'])\
				.merge(c_sample_tot_ints, on=['chr1', 'chr2'])\
				.assign(expected = lambda x: (x.coverage1*x.coverage2)/(x.total_coverage * 2))\
				.assign(oe = lambda x: x['count'] / x['expected'])\
				.drop(['coverage1', 'coverage2', 'total_coverage'], axis=1)
	return c_sample_compartment_agg_oe


def interchromosomal_oe(sample_compartment_agg):
	c_sample_compartment_agg = sample_compartment_agg[sample_compartment_agg.chr1 != sample_compartment_agg.chr2].copy()
	c_sample_comp1_coverage = c_sample_compartment_agg.groupby(['chr1', 'chr2', 'compartment1'], observed=True)['count']\
									.sum().to_frame('coverage1').reset_index()
	c_sample_comp2_coverage = c_sample_compartment_agg.groupby(['chr1', 'chr2', 'compartment2'], observed=True)['count']\
									.sum().to_frame('coverage2').reset_index()
	c_sample_tot_ints = c_sample_compartment_agg.groupby(['chr1', 'chr2'], observed=True)['count'].sum().to_frame('total_coverage').reset_index()
	c_sample_compartment_agg_oe = c_sample_compartment_agg\
					.merge(c_sample_comp1_coverage, on=['chr1', 'chr2', 'compartment1'])\
					.merge(c_sample_comp2_coverage, on=['chr1', 'chr2', 'compartment2'])\
					.merge(c_sample_tot_ints, on=['chr1', 'chr2'])\
					.assign(expected = lambda x: (x.coverage1*x.coverage2)/(x.total_coverage))\
					.assign(oe = lambda x: x['count'] / x['expected'])\
					.drop(['coverage1', 'coverage2', 'total_coverage'], axis=1)
	return c_sample_compartment_agg_oe



def load_sample_compartment_aggregations(sample_path):
	sample_compartment_agg = pd.read_csv(sample_path, sep='\t')
	sample_compartment_agg = sample_compartment_agg.rename(columns = {'category1': 'compartment1', 'category2': 'compartment2'})
	sample_compartment_agg['chr1'] = pd.Categorical(sample_compartment_agg['chr1'], categories=CHROMS, ordered=True)
	sample_compartment_agg['chr2'] = pd.Categorical(sample_compartment_agg['chr2'], categories=CHROMS, ordered=True)
	sample_compartment_agg['compartment1'] = pd.Categorical(sample_compartment_agg['compartment1'], categories=COMPARTMENTS_8_LABELS, ordered=True)
	sample_compartment_agg['compartment2'] = pd.Categorical(sample_compartment_agg['compartment2'], categories=COMPARTMENTS_8_LABELS, ordered=True)
	sample_compartment_agg_intra = intrachromosomal_oe(sample_compartment_agg)
	sample_compartment_agg_inter = interchromosomal_oe(sample_compartment_agg)
	sample_compartment_agg_oe = pd.concat([sample_compartment_agg_intra, sample_compartment_agg_inter], axis=0, ignore_index=True)
	return sample_compartment_agg, sample_compartment_agg_oe


def plot_comp_oe_intra_aggregate_heatmap(sample_compartment_agg_oe, output_path):
	X = sample_compartment_agg_oe[(sample_compartment_agg_oe.chr1 == sample_compartment_agg_oe.chr2) & \
								  (sample_compartment_agg_oe.compartment1 != 'Unknown') & \
								  (sample_compartment_agg_oe.compartment2 != 'Unknown')]\
						.groupby(['compartment1', 'compartment2'], observed=True)['oe'].mean()\
						.to_frame("avg_oe").reset_index()
	M = pd.pivot_table(X, index='compartment1', columns='compartment2', values='avg_oe')
	with sns.plotting_context('paper', font_scale=1.1):
		fig = plt.figure()
		sns.heatmap(M.T, cmap='bwr', vmin=0, vmax=2, 
					center=1, linewidth=0.1, square=True,
					annot=False, fmt=".1f", annot_kws={'fontsize':10},
					cbar_kws=dict(shrink=0.5, aspect=10, label='n. contacts / exp'))
		plt.yticks(rotation=0)
		fig.savefig(output_path, transparent=True, bbox_inches='tight')
		plt.close(fig)

def plot_comp_oe_inter_aggregate_heatmap(sample_compartment_agg_oe, output_path):
	X = sample_compartment_agg_oe[(sample_compartment_agg_oe.chr1 != sample_compartment_agg_oe.chr2) & \
								  (sample_compartment_agg_oe.compartment1 != 'Unknown') & \
								  (sample_compartment_agg_oe.compartment2 != 'Unknown')]\
						.groupby(['compartment1', 'compartment2'], observed=True)['oe'].mean()\
						.to_frame("avg_oe").reset_index()
	M = pd.pivot_table(X, index='compartment1', columns='compartment2', values='avg_oe')
	with sns.plotting_context('paper', font_scale=1.1):
		fig = plt.figure()
		sns.heatmap(M.T, cmap='bwr', vmin=0, vmax=2, 
					center=1, linewidth=0.1, square=True,
					annot=False, fmt=".1f", annot_kws={'fontsize':10},
					cbar_kws=dict(shrink=0.5, aspect=10, label='n. contacts / exp'))
		plt.yticks(rotation=0)
		fig.savefig(output_path, transparent=True, bbox_inches='tight')
		plt.close(fig)



def plot_comp_delta_intra_aggregate_heatmap(s12_delta, output_path):
	X = s12_delta[(s12_delta.chr1 == s12_delta.chr2) & \
				  (s12_delta.compartment1 != 'Unknown') & \
				  (s12_delta.compartment2 != 'Unknown')]
	Xother = X[X.compartment1 != X.compartment2].copy()			  
	comp1 = Xother.compartment1.copy()
	comp2 = Xother.compartment2.copy()

	Xother['compartment1'] = comp2
	Xother['compartment2'] = comp1

	X = pd.concat([X, Xother], axis=0)
	X = X.groupby(['compartment1', 'compartment2'], observed=True)['foldchange'].mean()\
						.to_frame("avg_oe").reset_index()
	M = pd.pivot_table(X, index='compartment1', columns='compartment2', values='avg_oe')
	with sns.plotting_context('paper', font_scale=1.1):
		fig = plt.figure()
		sns.heatmap(M.T.applymap(np.log2), cmap='RdBu_r', vmin=-1, vmax=1, 
					center=0, linewidth=0.1, square=True,
					annot=False, fmt=".1f", annot_kws={'fontsize':10},
					cbar_kws=dict(shrink=0.5, aspect=10, label="Fold-change"))
		plt.yticks(rotation=0)
		fig.savefig(output_path, transparent=True, bbox_inches='tight')
		plt.close(fig)
	return X


def plot_comp_delta_inter_aggregate_heatmap(s12_delta, output_path):
	X = s12_delta[(s12_delta.chr1 != s12_delta.chr2) & \
				  (s12_delta.compartment1 != 'Unknown') & \
				  (s12_delta.compartment2 != 'Unknown')].copy()

	Xother = X[X.chr1 != X.chr2].copy()			  
	chr1 = Xother.chr1.copy()
	comp1 = Xother.compartment1.copy()
	chr2 = Xother.chr2.copy()
	comp2 = Xother.compartment2.copy()

	Xother['chr1'] = chr2
	Xother['compartment1'] = comp2
	Xother['chr2'] = chr1
	Xother['compartment2'] = comp1

	X = pd.concat([X, Xother], axis=0)
	X = X.groupby(['compartment1', 'compartment2'], observed=True)['foldchange'].mean()\
						.to_frame("avg_oe").reset_index()
	M = pd.pivot_table(X, index='compartment1', columns='compartment2', values='avg_oe')
	with sns.plotting_context('paper', font_scale=1.1):
		fig = plt.figure()
		sns.heatmap(M.T.applymap(np.log2), cmap='RdBu_r', vmin=-0.3, vmax=0.3, 
					center=0, linewidth=0.1, square=True,
					annot=False, fmt=".1f", annot_kws={'fontsize':10},
					cbar_kws=dict(shrink=0.5, aspect=10, label="Fold-change"))
		plt.yticks(rotation=0)
		fig.savefig(output_path, transparent=True, bbox_inches='tight')
		plt.close(fig)
	return X


def main():
	args = process_arguments()

	os.makedirs(args.output_path, exist_ok=True)

	sample1_compartment_agg, sample1_compartment_agg_oe = load_sample_compartment_aggregations(args.sample1)
	sample1_compartment_agg_oe = sample1_compartment_agg_oe.assign(sample = args.name1)

	plot_comp_oe_intra_aggregate_heatmap(sample1_compartment_agg_oe, os.path.join(args.output_path, f"{args.name1}_OE_INTRA.pdf"))
	plot_comp_oe_inter_aggregate_heatmap(sample1_compartment_agg_oe, os.path.join(args.output_path, f"{args.name1}_OE_INTER.pdf"))

	sample2_compartment_agg, sample2_compartment_agg_oe = load_sample_compartment_aggregations(args.sample2)
	sample2_compartment_agg_oe = sample2_compartment_agg_oe.assign(sample = args.name2)

	plot_comp_oe_intra_aggregate_heatmap(sample2_compartment_agg_oe, os.path.join(args.output_path, f"{args.name2}_OE_INTRA.pdf"))
	plot_comp_oe_inter_aggregate_heatmap(sample2_compartment_agg_oe, os.path.join(args.output_path, f"{args.name2}_OE_INTER.pdf"))


	sample12_compartment_agg_oe = pd.concat([sample1_compartment_agg_oe, sample2_compartment_agg_oe], axis=0, ignore_index=True)
	all_compartment_aggs_tab = pd.pivot_table(sample12_compartment_agg_oe, index=['chr1', 'compartment1', 'chr2', 'compartment2'], columns='sample', values='oe')

	s12_delta = all_compartment_aggs_tab[args.name2] / all_compartment_aggs_tab[args.name1]
	s12_delta = s12_delta.to_frame('foldchange').reset_index()
	
	plot_comp_delta_intra_aggregate_heatmap(s12_delta, os.path.join(args.output_path, f"{args.name1}_VS_{args.name2}_FC_INTRA.pdf"))
	plot_comp_delta_inter_aggregate_heatmap(s12_delta, os.path.join(args.output_path, f"{args.name1}_VS_{args.name2}_FC_INTER.pdf"))


if __name__ == '__main__':
	main()



