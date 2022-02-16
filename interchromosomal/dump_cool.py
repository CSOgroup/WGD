#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Dumps Hi-C interactions at specified resolution and performs inter-chromosomal balancing.
The scripts dumps:
- Raw counts for each pairs of bins
- Balanced counts from IC correction for inter-chromosomal interactions
- Rescaled balanced counts from IC correction for inter-chromosomal interactions

This script takes as input:
- an Hi-C matrix path, in Cooler format
- the output path where to store the result, which is a 
 tsv file, having as columns [chr1, start1, end2, chr2, start2, end2, count, balanced, balanced_rescaled]
- [OPTIONAL] resolution: at which resolution to dump (defaults to 10Mb)

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
import logging
import os
import subprocess
from cooler import Cooler, balance_cooler
import pandas as pd
import numpy as np
from pybedtools.bedtool import BedTool
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

TMP_COARSENED_PATH = ".coarsened.cool"


HG19_CHROMS = ["chr" + str(x) for x in range(1, 23)] + ['chrX']
MM9_CHROMS = ["chr" + str(x) for x in range(1, 20)] + ['chrX']

CHROMS = None

def process_arguments():
	parser = argparse.ArgumentParser(description='Dump Hi-C data at specified resolutions')
	parser.add_argument("cool_path", type=str, help='Path to the cool file')
	parser.add_argument('output_path', type=str, help='Output path where to store the dumped result')
	parser.add_argument('--resolution', type=int, const=10000000, nargs='?', help='Minimum distance between bins to include in the calculation')
	parser.add_argument('--genome', type=str, const="hg19", nargs='?', help='Reference genome')
	args = parser.parse_args()
	return args

def set_logger():
	logger = logging.getLogger("Dump cool")
	handler = logging.StreamHandler()
	formatter = logging.Formatter('%(name)-12s %(message)s')
	handler.setFormatter(formatter)
	logger.addHandler(handler)
	logger.setLevel(logging.INFO)
	return logger

def pixel_to_matrix(cool_ints_inter, bins, value_col='balanced'):
	cool_ints_inter = cool_ints_inter.merge(bins.rename(columns = lambda x: x + '1'), on=['chr1', 'start1', 'end1'])\
						.merge(bins.rename(columns = lambda x: x + '2'), on=['chr2', 'start2', 'end2'])

	inter_m = np.zeros((bins.shape[0], bins.shape[0]))*np.NaN
	inter_m[cool_ints_inter.bin_id1.values, cool_ints_inter.bin_id2.values] = cool_ints_inter[value_col].values
	inter_m[cool_ints_inter.bin_id2.values, cool_ints_inter.bin_id1.values] = cool_ints_inter[value_col].values
	return inter_m

def rescale_margins(m, n_rescaling_iterations=10):
	mb = m.copy()
	for _ in range(n_rescaling_iterations):
		rm = np.nansum(mb, axis=1)[:,np.newaxis]
		mb = mb/rm
		cm = np.nansum(mb, axis=0)[np.newaxis, :]
		mb = mb/cm
	return mb


def main():
	logger = set_logger()

	args = process_arguments()
	input_cool = args.cool_path
	output_path = args.output_path
	resolution = args.resolution
	genome = args.genome

	logger.info(f"Input Hi-C: {input_cool}")
	logger.info(f"Resolution: {resolution}")
	logger.info(f"Output path: {output_path}")
	logger.info(f"Genome: {genome}")

	if genome == 'hg19':
		CHROMS = HG19_CHROMS
	elif genome == 'mm9':
		CHROMS = MM9_CHROMS
	else:
		raise ValueError("Unknown genome")



	cool = Cooler(input_cool)

	start_binsize = cool.binsize
	logger.info(f"Starting binsize: {start_binsize}")

	if start_binsize < resolution:
		factor = resolution//start_binsize
		logger.info(f"Coarsing the data by a factor of {factor}")
		subprocess.run(f"cooler coarsen -k {factor} -o {TMP_COARSENED_PATH} {input_cool}", shell=True)
		cool = Cooler(TMP_COARSENED_PATH)

	logger.info("Starting balancing")
	bins = cool.bins()[:].rename(columns={'chrom':'chr'})
	if not bins.iloc[0].chr.startswith("chr"):
		bins['chr'] = 'chr' + bins.chr.astype(str)
	chroms_to_remove = set(bins.chr).difference(CHROMS)
	blacklist = bins[bins.chr.isin(chroms_to_remove)].index.tolist()
	logger.info("Blacklisting {} bins from chroms {}".format(len(blacklist), ", ".join(chroms_to_remove)))

	logger.info("Balancing trans-only")
	_ = balance_cooler(cool, 
				   trans_only=True,
				   rescale_marginals=False,
				   ignore_diags=1,
				   store=True,
				   blacklist=blacklist,
				   store_name='IC_trans', 
				   min_nnz=10, 
				   tol=1e-09)

	logger.info("Loading result")
	cool_ints = cool.matrix(balance='IC_trans', as_pixels=True, join=True)[:]
	if not cool_ints.chrom1.iloc[0].startswith("chr"):
		cool_ints['chrom1'] = "chr" + cool_ints.chrom1.astype(str)
		cool_ints['chrom2'] = "chr" + cool_ints.chrom2.astype(str)
	cool_ints = cool_ints.rename(columns = {'chrom1': 'chr1', 'chrom2': 'chr2'})

	cool_ints = cool_ints[(cool_ints.chr1.isin(CHROMS)) & (cool_ints.chr2.isin(CHROMS))]
	bins = bins[bins.chr.isin(CHROMS)]
	bins['bin_id'] = np.arange(bins.shape[0], dtype = int)#bins.index

	logger.info(f"{bins}")

	logger.info("Rebalancing making marginals = 1")
	cool_ints_inter = cool_ints[cool_ints.chr1 != cool_ints.chr2]

	logger.info(cool_ints_inter)

	inter_m = pixel_to_matrix(cool_ints_inter, bins, value_col='balanced')
	inter_m_rescaled = rescale_margins(inter_m)
	inter_m_rescaled[np.tril_indices_from(inter_m_rescaled, -1)] = np.NaN

	bin_id1, bin_id2 = np.where(~np.isnan(inter_m_rescaled))
	values = inter_m_rescaled[bin_id1, bin_id2]

	cool_ints_inter_rescaled = pd.DataFrame({
		'bin_id1': bin_id1,
		'bin_id2': bin_id2,
		'balanced_rescaled': values
		})
	cool_ints_inter_rescaled = bins.rename(columns = lambda x: x + "1")\
									.merge(cool_ints_inter_rescaled, on=['bin_id1'])\
									.drop('bin_id1', axis=1)\
									.merge(bins.rename(columns = lambda x: x + "2"), on='bin_id2')\
									.drop('bin_id2', axis=1)
	cool_ints_inter = cool_ints_inter.merge(cool_ints_inter_rescaled, on=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])

	logger.info("Aggregating results together")
	cool_ints_final = pd.concat([
			cool_ints[cool_ints.chr1 == cool_ints.chr2].assign(balanced_rescaled = np.NaN, balanced=np.NaN),
			cool_ints_inter
		], axis=0, ignore_index=True)

	if os.path.isfile(TMP_COARSENED_PATH):
		logger.info(f"Removing {TMP_COARSENED_PATH}")
		os.remove(TMP_COARSENED_PATH)

	logger.info("Saving result")
	cool_ints_final.to_csv(output_path, sep='\t', index=False, header=True)

if __name__ == '__main__':
	main()


