import os
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from statsmodels.nonparametric.smoothers_lowess import lowess
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import seaborn as sns


def get_info(sample):
	if sample.startswith("RPE_TP53"):
		cell_line = "RPE_TP53"
		condition = sample.split("_")[2]
		target = sample.split("_")[3]
	elif sample.startswith("CPA-TP53"):
		cell_line = "CPA"
		condition = sample.split("_")[1]
		target = sample.split("_")[2]
	else:
		raise ValueError("Unknown situation")

	return cell_line, condition, target

def get_peaks_path(cell_line, target):
	for f in os.listdir(CHIPSEQ_PEAKS_PATH):
		if f.startswith(cell_line + "_" + target) and f.endswith("_consensusPeaks.tsv"):
			return os.path.join(CHIPSEQ_PEAKS_PATH, f)
	raise ValueError("Not found")



def plot_point_kde(ax, data, x, y, sample1="sample1", sample2="sample2", **kwargs):
	xx = data[x]
	yy = data[y]
	xy = np.vstack([xx,yy])
	z = gaussian_kde(xy)(xy)
	cax = ax.scatter(xx, yy, c=z, **kwargs)

	xmin, xmax = min(data[x].min(), data[y].min()), max(data[x].max(), data[y].max())
	plt.plot([xmin, xmax], [xmin, xmax], linestyle = '--', color = 'black', linewidth=2)
	z = lowess(data[y].values, data[x].values)
	plt.plot(z[:, 0], z[:, 1], color = 'red', linewidth=3)
	plt.xlabel(f"{sample1} signal\n(FC over input)")
	plt.ylabel(f"{sample2} signal\n(FC over input)")
	plt.title(f"Number of peaks = {data.shape[0]}")
	sns.despine(ax=ax)



OUTPUT_PATH = "figures/ChIPseq"

CHIPSEQ_PEAKS_PATH = "data/chipseq/peaks"

COMPARISONS = [
	("CPA-TP53-Clone3_Ctrl_CTCF vs CPA-TP53-Clone3_WGD_CTCF"),
	("CPA-TP53-Clone3_Ctrl_H3K9me3 vs CPA-TP53-Clone3_WGD_H3K9me3"),
	("RPE_TP53_Ctrl_CTCF vs RPE_TP53_WGD_CTCF"),
	("RPE_TP53_Ctrl_H3K9me3 vs RPE_TP53_WGD_H3K9me3"),
]

os.makedirs(OUTPUT_PATH, exist_ok=True)

print("Comparing ChIP-seq signals across conditions")
for comparison in COMPARISONS:
	sample1, sample2 = comparison.split(" vs ")
	name = os.path.join(OUTPUT_PATH, f"{sample1}_vs_{sample2}_signalComparison.pdf")
	if os.path.isfile(name):
		continue
	cell_line1, condition1, target1 = get_info(sample1)
	cell_line2, condition2, target2 = get_info(sample2)
	sample1_str = f"{sample1}"
	sample2_str = f"{sample2}"
	assert target1 == target2, "Different targets!"
	assert cell_line1 == cell_line2, "Different cell lines!"
	peaks_path = get_peaks_path(cell_line1, target1)
	print(f"{sample1_str:50}{sample2_str:50}{peaks_path}")
	peaks = pd.read_csv(peaks_path, sep="\t")


	sample1_signal = f"maxSignal_{sample1}"
	sample2_signal = f"maxSignal_{sample2}"


	fig, ax = plt.subplots(1, 1, figsize=(5, 5))
	plot_point_kde(ax, 
				   peaks, 
				   sample1_signal, 
				   sample2_signal, 
				   sample1=sample1,
				   sample2=sample2,
				   cmap='viridis', 
				   norm=LogNorm(),
				   rasterized=True)
	fig.savefig(name)
	plt.close(fig)
