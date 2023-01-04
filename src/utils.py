import os
import uuid
import subprocess
import pandas as pd


def load_chromsizes():
	CHROMS = ['chr' + str(x) for x in range(1, 23)] + ['chrX']
	CHROMSIZES = pd.read_csv("resources/hg19.main.chrom.sizes", sep='\t', header=None, names=['chr', 'length'])
	CHROMSIZES = CHROMSIZES[CHROMSIZES.chr.map(lambda x: x).isin(CHROMS)]
	return CHROMSIZES



def run_CBS(df, 
			smooth=True, 
			threshold=0.01, 
			nperm=10000, 
			perm_method = "perm", 
			undo_split = 'none',
			min_width=2):
	""" Runs the Circular Binary Segmentation algorithm.
	It accepts:
	- A dataframe having the columns: 
		chr: chromosome
		start: start of the bin
		ratio: value to segment. It is assumed a log2 ratio
		[weight]: weight of each bin. May be used to smooth 
				  the result (proportional to the bin variance)
		[...]: any other column is not considered
	- smooth: if True, it smoothes the ratio values before 
			  running CBS [true]
	- threshold: significance value to use to call the segments [0.01]
	- nperm: number of permutations for p-value calculation
	- perm_method: type of permutation method to calculate p-values
	- min_width: minimum number of bins for a changed segment

	It returns:
	- A dataframe having the columns 
		chr: chromosome
		start: starting position of the segment along the chromosome
		end: end position of the segment along the chromosome
		mean_value: average value of *ratio* in the segment
		n_points: number of bins which compose the segment 
	"""

	rand_input_file_name = str(uuid.uuid4())
	rand_output_file_name = str(uuid.uuid4())
	
	CBS_SCRIPT="src/cbs.R"
	
	CBS_COMMAND=f"Rscript {CBS_SCRIPT} {rand_input_file_name} {rand_output_file_name}"
	if smooth:
		CBS_COMMAND += " --smooth"
	CBS_COMMAND += f" --threshold {threshold}"
	CBS_COMMAND += f" --nperm {nperm}"
	CBS_COMMAND += f" --perm_method {perm_method}"
	CBS_COMMAND += f" --undo_split {undo_split}"
	CBS_COMMAND += f" --min_width {min_width}"	
	# print(CBS_COMMAND)
	# CBS_COMMAND += " 2> /dev/null"

	df.to_csv(rand_input_file_name, sep='\t', index=False, header=True)
	subprocess.run(CBS_COMMAND, shell=True)
	result = pd.read_csv(rand_output_file_name, sep="\t")
	os.remove(rand_input_file_name)
	os.remove(rand_output_file_name)
	result = result.rename(columns={'chrom': 'chr', 
									'loc.start': 'start', 
									'loc.end': 'end', 
									'seg.mean': 'mean_value', 
									'num.mark': 'n_points'})
	result = result[['chr', 'start', 'end', 'mean_value', 'n_points']]
	return result