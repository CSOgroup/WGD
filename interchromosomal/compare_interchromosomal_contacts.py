#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Compares the inter-chromosomal contacts between two Hi-C samples.
It gives as output a set of plots showing the differences in relative interactions
between pairs of chromosomes between the two samples.

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
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation


CHROMS = ['chr' + str(x) for x in range(1, 23)] + ['chrX']
CHROMSIZES = pd.read_csv("resources/hg19.main.chrom.sizes", sep='\t', header=None, names=['chr', 'length'])
CHROMSIZES = CHROMSIZES[CHROMSIZES.chr.map(lambda x: x).isin(CHROMS)]
CHROMSIZES = CHROMSIZES.sort_values('length', ascending=False)
CHROMS_BY_SIZE = CHROMSIZES.chr.tolist()

LONG_CHROMS =  CHROMS[:14] + [CHROMS[-1]]
SHORT_CHROMS = CHROMS[14:-1]


def process_arguments():
    parser = argparse.ArgumentParser(description='Compares the inter-chromosomal contacts between two Hi-C experiments')
    parser.add_argument("sample1", type=str, help='Path to the first sample (control)')
    parser.add_argument('sample2', type=str, help='Path to the second sample (treatment)')
    parser.add_argument("output_path", type=str, help="Where to store the results")
    parser.add_argument("--name1", type=str, default="sample1", help="Name of first sample")
    parser.add_argument("--name2", type=str, default="sample2", help="Name of second sample")
    args = parser.parse_args()
    return args


def load_hic_dumps(spath, chrom_order = CHROMS_BY_SIZE):
    s_ints = pd.read_csv(spath, sep="\t")
    s_ints['bin_id1'] = s_ints.chr1 + ":" + s_ints.start1.astype(str) + "-" + s_ints.end1.astype(str)
    s_ints['bin_id2'] = s_ints.chr2 + ":" + s_ints.start2.astype(str) + "-" + s_ints.end2.astype(str)
    s_ints_other = s_ints[s_ints.bin_id1 != s_ints.bin_id2].copy()
    x1 = s_ints_other[['chr1', 'start1', 'end1', 'bin_id1']].copy()
    x2 = s_ints_other[['chr2', 'start2', 'end2', 'bin_id2']].copy()
    s_ints_other.loc[:, ['chr1', 'start1', 'end1', 'bin_id1']] = x2.values
    s_ints_other.loc[:, ['chr2', 'start2', 'end2', 'bin_id2']] = x1.values
    s_ints_both = pd.concat([s_ints, s_ints_other], axis=0, ignore_index=True)
    s_ints_both['chr1'] = pd.Categorical(s_ints_both['chr1'], categories=chrom_order, ordered=True)
    s_ints_both['chr2'] = pd.Categorical(s_ints_both['chr2'], categories=chrom_order, ordered=True)
    s_ints_both['chr1_type'] = s_ints_both.chr1.map(lambda x: 'L' if x in LONG_CHROMS else 'S')
    s_ints_both['chr2_type'] = s_ints_both.chr2.map(lambda x: 'L' if x in LONG_CHROMS else 'S')
    s_ints_both['chr_type_int'] = s_ints_both.apply(lambda x: "".join(sorted([x.chr1_type, x.chr2_type])), axis=1)
    return s_ints_both



def plot_inter_chromosomal_map_chrom_level_oe(g, output_path, chroms=CHROMS_BY_SIZE):
    M = pd.pivot_table(g, 
                       index='chr1', 
                       columns='chr2', 
                       values='balanced_rescaled', 
                       observed=True)
    M = M.reindex(index = chroms, columns = chroms)
    
    with sns.plotting_context('paper', font_scale=1.2):
        fig = plt.figure()
        sns.heatmap(M, cmap='bwr', square=True,
                    center = 1/23, vmin=0.5/23, vmax=1.5/23,
                    linewidth=0.1,
                    cbar_kws={'label': 'balanced counts', 'shrink': 0.5, 'aspect': 8})
        fig.savefig(output_path, transparent=True, bbox_inches='tight')
        plt.close(fig)



def plot_inter_chromosomal_map_chrom_level_foldchange(delta, output_path, chroms=CHROMS_BY_SIZE):
    M = pd.pivot_table(delta[delta.chr1 != delta.chr2], 
            index='chr1', 
            columns='chr2', 
            values='foldchange', 
            observed=True)
    M = M.reindex(index = chroms, columns = chroms)
    with sns.plotting_context('paper', font_scale=1.2):
        fig = plt.figure()
        sns.heatmap(M.applymap(lambda x: np.log2(x)), 
                    vmin=-0.5, vmax=0.5, 
                    cmap='RdBu_r',
                    linewidth=0.1,
                    cbar_kws={'label': '$log_2$(fold-change)\nof balanced counts', 'shrink': 0.5, 'aspect': 8}, 
                    square=True)
        fig.savefig(output_path, transparent=True, bbox_inches='tight')
        plt.close(fig)


def plot_inter_chromosomal_chrom_level_raincloud(delta, output_path):
    with sns.plotting_context('paper', font_scale = 1.2):
        fig, ax = plt.subplots(1, 1, figsize=(3, 3))

        sns.violinplot(x = 'chr_type_int', y = 'foldchange', data = delta, palette = 'colorblind',
                       order = ['LL', 'LS', 'SS'], ax = ax)
        sns.stripplot(ax = ax, data=delta, x = 'chr_type_int', y = 'foldchange', color='black',
                    order = ['LL', 'LS', 'SS'], size=2)
        ax.set_ylim(0.6, 1.45)
        add_stat_annotation(ax, data=delta, x='chr_type_int', y='foldchange', order = ['LL', 'LS', 'SS'],
                            box_pairs = [('LL', 'LS'),
                                         ('LL', 'SS'),
                                         ('LS', 'SS')], test='Mann-Whitney', loc='outside',
                                         text_format = 'full')
        sns.despine()
        plt.xlabel("Chromosome interactions")
        plt.ylabel('$log_2$(fold-change)\nof balanced counts')
        fig.savefig(output_path, transparent=True, bbox_inches='tight')
        plt.close(fig)



def main():
    args = process_arguments()

    os.makedirs(args.output_path, exist_ok=True)

    sample1_both = load_hic_dumps(args.sample1, chrom_order = CHROMS).assign(sample = args.name1)
    sample2_both = load_hic_dumps(args.sample2, chrom_order = CHROMS).assign(sample = args.name2)

    plot_inter_chromosomal_map_chrom_level_oe(sample1_both, os.path.join(args.output_path, f"{args.name1}_OE.pdf"))
    plot_inter_chromosomal_map_chrom_level_oe(sample2_both, os.path.join(args.output_path, f"{args.name2}_OE.pdf"))

    sample12 = pd.concat([sample1_both, sample2_both], axis=0, ignore_index=True)

    all_ints_tab = pd.pivot_table(sample12, 
                      index=['chr1', 'start1', 'end1', 'chr1_type', 'chr2', 'start2', 'end2', 'chr2_type', 'chr_type_int'], 
                      columns = 'sample', 
                      values='balanced_rescaled', 
                      observed=True)

    delta = all_ints_tab[args.name2] / all_ints_tab[args.name1]
    delta = delta.to_frame('foldchange').reset_index()

    plot_inter_chromosomal_map_chrom_level_foldchange(delta, os.path.join(args.output_path, f"{args.name1}_VS_{args.name2}_FC.pdf"))
    plot_inter_chromosomal_chrom_level_raincloud(delta, os.path.join(args.output_path, f"{args.name1}_VS_{args.name2}_RainCloud.pdf"))


if __name__ == '__main__':
    main()




