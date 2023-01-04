import os
import pandas as pd
import numpy as np
from cooler import Cooler
from cooler.fileops import list_coolers
from cooltools.lib.numutils import observed_over_expected
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from pybedtools.bedtool import BedTool
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append('./src/')
from utils import run_CBS, load_chromsizes

CHROMS = [f'chr{x}' for x in range(1, 23)] + ['chrX']
CHROMSIZES = load_chromsizes()
CHROMS_BY_SIZE = CHROMSIZES.sort_values('length', ascending=False).chr.tolist()

BINSIZE=1000000
COVERAGE_BINSIZE=5000000
EXCLUDED_CHROMS=["chrY", 'chrM']
MAX_DEVIATION_FROM_GENOMEWIDE_MEDIAN = 0.4

GAPS = pd.read_csv("resources/hg19_gap.txt", sep="\t", header=None, usecols=[0, 1, 2, 6], names = ['chrom', 'start', 'end', 'gap_type'])
GAPS = GAPS[GAPS.gap_type.isin(['telomere', 'centromere'])]

UNSEGREGATED_CELLS = [
 'RPE_WGD_schic_003',
 'RPE_WGD_schic_005',
 'RPE_WGD_schic_014',
 'RPE_WGD_schic_017',
 'RPE_WGD_schic_013',
 'RPE_WGD_schic_018',
 'RPE_WGD_schic_008',
 'RPE_WGD_schic_023'
]

CELL_CATEGORY_PALETTE = {
	"Control": "grey",
	"WGD without LCS": "black",
	"WGD with LCS": "red"
}


CHROMOSOME_INTERACTIONS_PALETTE = {
	"LL": "blue",
	"LS": "red",
	"SS": "lightblue"
}

def to_cell_category(x):
	if x.condition == "RPE_TP53":
		return "Control"
	else:
		if x.cell in UNSEGREGATED_CELLS:
			return "WGD with LCS"
		else:
			return "WGD without LCS"

def make_symmetric(all_chrom_interactions):
	x = all_chrom_interactions[all_chrom_interactions.chr1 != all_chrom_interactions.chr2].copy()
	x1 = x.chr1.copy()
	x2 = x.chr2.copy()
	x['chr1'] = x2
	x['chr2'] = x1 
	all_chrom_interactions_both = pd.concat([
			all_chrom_interactions,
			x,
		], axis=0, ignore_index=True)
	return all_chrom_interactions_both

def get_AB_compartments(bulk_compartments_path):
	bulk_compartments = pd.read_csv(bulk_compartments_path, sep="\t", header=None, 
									names = ['chr', 'start', 'end', 'compartment_label', 
											 "compartment_label_8", "name", "tickStart", 
											 "tickEnd", "color"])
	bulk_compartments['compartment_label_2'] = bulk_compartments.compartment_label.map(lambda x: x.split(".")[0])
	bulk_compartments = bulk_compartments.sort_values(['chr', 'start', 'end']).reset_index(drop=True)

	bulk_compartments['consecutive'] = ((bulk_compartments.chr != bulk_compartments.chr.shift()) | \
						   			 (bulk_compartments.start > bulk_compartments.end.shift() + 1000000) | 
						   			 (bulk_compartments["compartment_label_2"] != bulk_compartments["compartment_label_2"].shift())).cumsum()

	bulk_compartments = bulk_compartments.groupby('consecutive', as_index=False, sort=False)\
				     .apply(lambda x: pd.Series({'chr': x.chr.iloc[0], 
				   							     'start': x.start.min(), 
				   							     'end': x.end.max(), 
				   							     'comp': x["compartment_label_2"].iloc[0]}))\
				     .drop("consecutive", axis=1)\
				     .assign(comp = lambda x: pd.Categorical(x.comp, categories=["A", "B"], ordered=True),
				     		 chr = lambda x: pd.Categorical(x.chr, categories=CHROMS, ordered=True))
	bulk_compartments = bulk_compartments.sort_values(['chr', 'start', 'end']).reset_index(drop=True)
	bulk_compartments['domain_id'] = bulk_compartments.index
	return bulk_compartments


def call_schic_compartments(cell_path, chrom, bulk_compartments, phasing_track_overlap=0.6, n_comps=10):
	H_cool = Cooler(cell_path)
	if chrom not in H_cool.chromnames:
		xc = chrom[3:]
	else:
		xc = chrom

	bins = H_cool.bins().fetch(xc)
	if not bins.iloc[0].chrom.startswith('chr'):
		bins['chrom'] = "chr" + bins.chrom.astype(str)
	else:
		bins['chrom'] = bins.chrom.astype(str)
	
	bins = bins.rename(columns = {'chrom': 'chr'})
	bins['bin_id'] = bins.index
	bins = bins[['chr', 'start', 'end', 'bin_id']]
	bins = BedTool.from_dataframe(bins)\
				.sort()\
				.map(BedTool.from_dataframe(bulk_compartments).sort(), 
					 c=4, o='distinct', 
					 null='Unknown', 
					 f=phasing_track_overlap)\
				.to_dataframe(names = bins.columns.tolist() + ['comp'])
	H = H_cool.matrix(balance=False).fetch(xc)

	n_interactions = np.triu(H).sum()
	n_single_interactions = np.count_nonzero(np.triu(H))

	H[~np.isfinite(H)] = 0
	mask = H.sum(axis=0) > 0
	OE, _, _, _ = observed_over_expected(H, mask)
	OE -= 1.0
	OE[~mask, :] = 0
	OE[:, ~mask] = 0
	OE_corr = pd.DataFrame(OE, index=bins.bin_id, columns=bins.bin_id).corr(method='pearson')

	bad_bins = np.isnan(OE_corr).all(axis=1)
	bins = bins.merge(bad_bins.to_frame('is_bad'), on='bin_id')

	good_bins = bins[~bins.is_bad]
	OE_corr = OE_corr.reindex(index=good_bins.bin_id.values, columns=good_bins.bin_id.values)

	pca = PCA(n_components = n_comps) 
	pca.fit(OE_corr) 
	X_t = pca.transform(OE_corr)
	X_t = pd.DataFrame(X_t, columns = [f"comp{x}" for x in range(n_comps)], index=good_bins.bin_id)
	X_t = good_bins.merge(X_t.reset_index(), on='bin_id')
	return X_t


def get_cell_interactions(cell_cool, excluded_chromosomes=None):
	excluded_chromosomes = [] if excluded_chromosomes is None else excluded_chromosomes
	cell_ints = cell_cool.matrix(balance=None, as_pixels=True, join=True)[:]
	cell_ints_other = cell_ints.copy()
	x1 = cell_ints[['chrom1', 'start1', 'end1']].copy()
	x2 = cell_ints[['chrom2', 'start2', 'end2']].copy()
	cell_ints_other.loc[:, ['chrom1', 'start1', 'end1']] = x2.values
	cell_ints_other.loc[:, ['chrom2', 'start2', 'end2']] = x1.values
	cell_ints_both = pd.concat([cell_ints, cell_ints_other], axis=0, ignore_index=True).drop_duplicates()
	cell_ints_both = cell_ints_both.sort_values(['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']).reset_index(drop=True)
	cell_ints_both = cell_ints_both[(~cell_ints_both.chrom1.isin(excluded_chromosomes)) & (~cell_ints_both.chrom2.isin(excluded_chromosomes))]
	return cell_ints_both


def get_hic_coverage(cell_cool, binsize=None, excluded_chromosomes=None):
	excluded_chromosomes = [] if excluded_chromosomes is None else excluded_chromosomes
	cell_ints_both = get_cell_interactions(cell_cool, excluded_chromosomes=excluded_chromosomes)
	cell_bins = cell_cool.bins()[:]
	cell_bins = cell_bins[(~cell_bins.chrom.isin(excluded_chromosomes))]
	coverage = cell_bins.merge(
						cell_ints_both.groupby(['chrom1', 'start1', 'end1'])['count']\
							 .sum()\
							 .to_frame("coverage")\
							 .reset_index()\
							 .rename(columns = {'chrom1': 'chrom', 'start1': 'start', 'end1': 'end'}),
					on = ['chrom', 'start', 'end'], how = "left")\
					.assign(coverage = lambda x: x.coverage.fillna(0).astype(int))
	if binsize is not None:
		coverage = coverage.assign(start = lambda x: (x.start//binsize)*binsize)\
			 .groupby(["chrom", 'start'])["coverage"]\
			 .sum()\
			 .to_frame("coverage")\
			 .reset_index()\
			 .assign(end = lambda x: x.start + binsize)\
			 [['chrom', 'start', 'end', 'coverage']]

	median_coverage = coverage.coverage.median()
	coverage['genomewideNormalized_coverage'] = coverage['coverage'] / median_coverage

	chromosome_coverage = coverage.groupby("chrom")['coverage'].median().to_frame("chromosome_coverage").reset_index()
	coverage = coverage.merge(chromosome_coverage, on = "chrom").assign(chromosomeNormalized_coverage = lambda x: x['coverage'] / x['chromosome_coverage'])
	coverage = coverage.drop("chromosome_coverage", axis=1)
	return coverage

def call_cnvs(cell_coverage, binsize, ratio_col = "chromosomeNormalized_coverage", gaps = None):
	filtered_coverage = cell_coverage[cell_coverage[ratio_col] != 0]

	if gaps is not None:
		filtered_coverage = BedTool.from_dataframe(filtered_coverage)\
								   .sort()\
								   .subtract(BedTool.from_dataframe(gaps).sort(), A=True)\
								   .to_dataframe(names = filtered_coverage.columns)

	segments = run_CBS(filtered_coverage.assign(ratio = lambda x: x[ratio_col].map(np.log2))[["chrom", "start", 'end', "ratio"]], min_width=2)
	segments = segments.rename(columns = {'chr': 'chrom'})
	segments['end'] += binsize
	segments = BedTool.from_dataframe(segments)\
					.sort()\
					.map(BedTool.from_dataframe(cell_coverage).sort(), 
						 c = (np.where(~cell_coverage.columns.isin(['chrom', 'start', 'end']))[0] + 1).tolist(), 
						 o = 'median')\
					.to_dataframe(names = segments.columns.tolist() + [c for c in cell_coverage.columns if c not in segments.columns])
	return segments



def plot_cell_segregation_scores(inter_scores):
	sns.stripplot(data=inter_scores, x='condition', y='inter_score', hue = "cell_category", palette = CELL_CATEGORY_PALETTE, ax=ax, zorder=1)
	sns.boxplot(data=inter_scores, x='condition', y='inter_score', 
				showfliers=False, boxprops = {'facecolor':'none', 'edgecolor':'black'}, 
				ax=ax, zorder=0)
	ax.set_ylim(0.5, 1)
	ax.set_ylabel("Single-cell Inter-chr LCS [LS / (LL + SS)]")
	ax.set_xlabel("")
	ax.legend(bbox_to_anchor=(1, 1))
	sns.despine(ax = ax, trim=True)

def plot_inter_chromosomal_rank_difference(X):
	g = sns.jointplot(data = X, x = "Control", y="WGD with LCS", hue = "chr_type_int", palette = CHROMOSOME_INTERACTIONS_PALETTE,
					  xlim=(30, 250), ylim=(30, 250), marginal_kws = dict(common_norm=False, hue_order = ['LS', "SS"]), marginal_ticks=False)
	g.set_axis_labels("Inter-chr interaction rank (Control)", "Inter-chr interaction rank (LCS-WGD)")
	g.ax_marg_x.axes.xaxis.set_visible(False)
	g.ax_marg_y.axes.xaxis.set_visible(False)
	return g

def plot_corr_heatmap(ax, M_control):
	sns.heatmap(M_control, 
				cmap='bwr', 
				vmin=-1, 
				vmax=1,
				linewidth=0.1,
				cbar_kws={'label': 'Pearson correlation of O/E', 'shrink': 0.5, 'aspect': 8}, 
				square=True,
				ax = ax)
	ax.set_xlabel("")
	ax.set_ylabel("")

def plot_corr_foldchange_heatmap(ax, M_wgd, M_control):
	sns.heatmap((M_wgd + 1) / (M_control + 1), 
				cmap='BrBG_r', 
				vmin=0, 
				vmax=2,
				linewidth=0.1,
				cbar_kws={'label': '$(corr_{WGD} + 1) / (corr_{Control} + 1)$', 'shrink': 0.5, 'aspect': 8}, 
				square=True,
				ax = ax)
	ax.set_xlabel("")
	ax.set_ylabel("")

def plot_cell_compartment_scores(all_cell_compartment_scores):
	g = sns.displot(data = all_cell_compartment_scores, 
				x = "silhouette_score_comp", 
				hue = "cell_category", 
				palette = CELL_CATEGORY_PALETTE,
				kind = 'kde', 
				common_norm=False)
	plt.xlabel("Single-cell compartment segregation")
	plt.ylabel("Density")
	return g

def plot_cnv_heatmap(heatmap, all_cells_cnv_covered_genome, limits):
	g = sns.clustermap(heatmap.applymap(lambda x: 1 if abs(x - 1) < MAX_DEVIATION_FROM_GENOMEWIDE_MEDIAN else x), 
						cmap='bwr', 
						vmin=0, 
						vmax=2, 
						xticklabels=False, 
						yticklabels=False, 
						row_cluster=False, 
						col_cluster=False,
						row_colors = pd.Series(palette, name = "Cell\npopulation"), 
						figsize=(30, 10),
						cbar_kws=dict(label = "Coverage ratio"),
						dendrogram_ratio = (0.1, 0))
	y = all_cells_cnv_covered_genome.set_index("cell")['cnv_perc_genome'].reindex(index = cell_order).values
	x = np.arange(len(y))
	g.ax_row_dendrogram.barh(x, y, color = 'black', align="edge", height=0.9)
	g.ax_row_dendrogram.set_ylim(0, 58)
	g.ax_row_dendrogram.invert_yaxis()
	g.ax_row_dendrogram.invert_xaxis()
	g.ax_row_dendrogram.text(0.5, -.02, "Fraction of altered genome", transform=g.ax_row_dendrogram.transAxes, horizontalalignment='center')
	for l in limits['min'].values.tolist() + [limits['max'].max()]:
		g.ax_heatmap.axvline(l, color = 'black')

	for t in limits.itertuples():
		g.ax_heatmap.text(t.mean, -1, t.Index, horizontalalignment='center', rotation=45, size=15)

	g.cax.set_visible(True)
	g.cax.set_aspect(3)
	g.ax_col_dendrogram.set_visible(False)


	g.ax_heatmap.set_xlabel("Genomic position")
	g.ax_heatmap.set_ylabel("Cell")
	return g


OUTPUT_PATH = "figures/ScHiC"
os.makedirs(OUTPUT_PATH, exist_ok=True)

SCHIC_PATH = "data/schic"


print("Analysing Sc-HiC inter-chromosomal contacts")
schic_samples = pd.read_excel(os.path.join(SCHIC_PATH, "schic_samples.xlsx"))
good_cells = schic_samples[schic_samples.flag.isnull()].sample_name.tolist()

RPE_chrom_chrom_interactions = pd.read_csv(os.path.join(SCHIC_PATH, "inter_chromosomal", "RPE_chromosome_interactions.tsv"), sep="\t")
RPE_chrom_chrom_interactions['chr1_type'] = RPE_chrom_chrom_interactions.chr1_type.map({"long": "L", "short": "S"})
RPE_chrom_chrom_interactions['chr2_type'] = RPE_chrom_chrom_interactions.chr2_type.map({"long": "L", "short": "S"})
RPE_chrom_chrom_interactions['chr_type_int'] = RPE_chrom_chrom_interactions.apply(lambda x: "".join(sorted([x.chr1_type, x.chr2_type])), axis=1)

metric = 'ice_count'

inter_scores = RPE_chrom_chrom_interactions[(RPE_chrom_chrom_interactions.chr1 != RPE_chrom_chrom_interactions.chr2)]\
	.groupby(['cell', 'condition', 'chr_type_int'])[metric]\
	.agg(lambda x: np.sum(x))\
	.to_frame("count")\
	.reset_index()\
	.pivot(index=['cell', 'condition'], columns='chr_type_int', values='count')\
	.assign(inter_score = lambda x: x['LS'] / (x['LL'] + x['SS']) )\
	.reset_index()
inter_scores['cell_category'] = inter_scores.apply(to_cell_category, axis=1)


fig, ax = plt.subplots(1, 1, figsize=(3, 5))
plot_cell_segregation_scores(inter_scores[inter_scores.cell.isin(good_cells)])
fig.savefig(os.path.join(OUTPUT_PATH, "inter_chromosomal_segregation.pdf"))
plt.close(fig)



print("Analysing inter-chromosomal interacations rankings between Control and WGD with LCS cells")
RPE_chrom_chrom_interactions = RPE_chrom_chrom_interactions[RPE_chrom_chrom_interactions.cell.isin(good_cells)].copy()

RPE_chrom_chrom_interactions["rank"] = RPE_chrom_chrom_interactions\
											.groupby(["cell"])\
											['ice_count']\
											.rank()
RPE_chrom_chrom_interactions['cell_category'] = RPE_chrom_chrom_interactions.apply(to_cell_category, axis=1)

X = RPE_chrom_chrom_interactions\
			.groupby(["cell_category", "chr1", 'chr2', 'chr_type_int'])\
			['rank']\
			.mean()\
			.dropna()\
			.to_frame("score")\
			.reset_index()\
			.pivot(index = ['chr1', 'chr2', 'chr_type_int'], columns = 'cell_category', values = 'score')\
			.reset_index()

g = plot_inter_chromosomal_rank_difference(X)
g.savefig(os.path.join(OUTPUT_PATH, "inter_chromosomal_rankDifference.pdf"))
plt.close(g.fig)

print("Pseudo-bulk inter-chromosomal analysis of ScHi-C data")

pseudo_chrom_interactions = RPE_chrom_chrom_interactions.groupby(["condition", "chr1", "chr2", "chr_type_int"])['raw_count'].sum().to_frame("count").reset_index()
pseudo_chrom_interactions_both = make_symmetric(pseudo_chrom_interactions)

pseudo_chr_coverage = pseudo_chrom_interactions_both[pseudo_chrom_interactions_both.chr1 != pseudo_chrom_interactions_both.chr2].groupby(['condition', 'chr1'])['count'].sum().reset_index().rename(columns = {'chr1': 'chr'})
pseudo_tot_coverage = pseudo_chrom_interactions[pseudo_chrom_interactions.chr1 != pseudo_chrom_interactions.chr2].groupby('condition')['count'].sum().reset_index()

pseudo_chrom_interactions_both = pseudo_chrom_interactions_both.merge(pseudo_chr_coverage.rename(columns = {'chr': 'chr1', "count": "coverage1"}), on = ["condition", "chr1"] )\
	.merge(pseudo_chr_coverage.rename(columns = {'chr': 'chr2', "count": "coverage2"}), on = ["condition", "chr2"] )\
	.merge(pseudo_tot_coverage.rename(columns = {'count': 'tot_coverage'}), on='condition')
pseudo_chrom_interactions_both = pseudo_chrom_interactions_both[pseudo_chrom_interactions_both.chr1 != pseudo_chrom_interactions_both.chr2]

pseudo_chrom_interactions_both['expected_count'] = (pseudo_chrom_interactions_both['coverage1']*pseudo_chrom_interactions_both['coverage2'])/(2*pseudo_chrom_interactions_both['tot_coverage'])
pseudo_chrom_interactions_both['oe_count'] = pseudo_chrom_interactions_both['count'] / pseudo_chrom_interactions_both['expected_count']


M_control = pd.pivot_table(pseudo_chrom_interactions_both[pseudo_chrom_interactions_both.condition == 'RPE_TP53'], index='chr1', columns = 'chr2', values='oe_count')
M_control = M_control.reindex(index=CHROMS_BY_SIZE, columns=CHROMS_BY_SIZE)
M_control = M_control.corr()

M_wgd = pd.pivot_table(pseudo_chrom_interactions_both[pseudo_chrom_interactions_both.condition == 'RPE_WGD'], index='chr1', columns = 'chr2', values='oe_count')
M_wgd = M_wgd.reindex(index=CHROMS_BY_SIZE, columns=CHROMS_BY_SIZE)
M_wgd = M_wgd.corr()

fig, ax = plt.subplots(1, 1)
plot_corr_heatmap(ax, M_control)
fig.savefig(os.path.join(OUTPUT_PATH, "RPE_TP53_pseudoBulk_Pearson.pdf"), transparent=True, bbox_inches='tight')
plt.close(fig)

fig, ax = plt.subplots(1, 1)
plot_corr_heatmap(ax, M_wgd)
fig.savefig(os.path.join(OUTPUT_PATH, "RPE_WGD_pseudoBulk_Pearson.pdf"), transparent=True, bbox_inches='tight')
plt.close(fig)

fig, ax = plt.subplots(1, 1)
plot_corr_foldchange_heatmap(ax, M_wgd, M_control)
fig.savefig(os.path.join(OUTPUT_PATH, "RPE_pseudoBulk_Pearson_foldChange.pdf"), transparent=True, bbox_inches='tight')
plt.close(fig)

print("Analysing single-cell Hi-C compartments")
os.makedirs(os.path.join(SCHIC_PATH, "compartments"), exist_ok=True)

bulk_compartments_path = "data/hic/compartment_domains/RPE_TP53_Ctrl_mega_CompartmentDomains.bed"
scool_path = "data/schic/maps/schic_matrices_binsize_1000000.scool"
bulk_compartments = get_AB_compartments(bulk_compartments_path)

phasing_track_overlap=0.6
n_comps = 2

all_cell_compartments_path = os.path.join(SCHIC_PATH, "compartments", "all_cell_compartments.tsv")
if not os.path.isfile(all_cell_compartments_path):
	all_cell_compartments = []
	for cell in list_coolers(scool_path):	
		cell_name = cell.split("/")[-1].split(".")[0].split("_binsize_")[0]
		print(f"Cell: {cell_name}")
		cell_path = scool_path + "::" + cell

		if not pd.isnull(schic_samples.set_index("sample_name").loc[cell_name, 'flag']):
			continue

		cell_compartments = []
		for chrom in CHROMS:
			print(f"\tChromosome: {chrom}")
			chrom_comps = call_schic_compartments(cell_path, chrom, bulk_compartments, phasing_track_overlap=phasing_track_overlap, n_comps=n_comps)
			cell_compartments.append(chrom_comps)
		cell_compartments = pd.concat(cell_compartments, axis=0, ignore_index=True)
		cell_compartments['cell'] = cell_name
		all_cell_compartments.append(cell_compartments)
	all_cell_compartments = pd.concat(all_cell_compartments, axis=0, ignore_index=True)
	all_cell_compartments.to_csv(all_cell_compartments_path, sep="\t", index=False, header=True)
else:
	all_cell_compartments = pd.read_csv(all_cell_compartments_path, sep="\t")

all_cell_compartment_scores = []
all_conversions = []
for gc, df in all_cell_compartments.groupby(['cell', 'chr']):
	cell, chrom = gc
	print(f"{cell} - {chrom}")
	X = df[df.comp != 'Unknown']
	silhouette_comp = silhouette_score(X[['comp0', 'comp1']], X['comp'], metric='euclidean')	
	all_cell_compartment_scores.append({'cell': cell, 
										'chr': chrom, 
										'silhouette_score_comp': silhouette_comp})
all_cell_compartment_scores = pd.DataFrame.from_dict(all_cell_compartment_scores)
all_cell_compartment_scores = all_cell_compartment_scores.merge(schic_samples.rename(columns = {'sample_name': 'cell'}), on='cell')
all_cell_compartment_scores['condition'] = all_cell_compartment_scores['cell_line'] + "_" + all_cell_compartment_scores['condition']
all_cell_compartment_scores['cell_category'] = all_cell_compartment_scores.apply(to_cell_category, axis=1)

g = plot_cell_compartment_scores(all_cell_compartment_scores)
g.savefig(os.path.join(OUTPUT_PATH, "singleCell_compartmentSegregation.pdf"))
plt.close(g.fig)

print("Analysing CNVs on ScHI-C data")
os.makedirs(os.path.join(SCHIC_PATH, "cnvs"), exist_ok=True)

name = os.path.join(SCHIC_PATH, "cnvs", "schic_cell_cnvs.tsv")
if not os.path.isfile(name):
	all_cell_cnvs = []
	for cell in list_coolers(scool_path):
		cell_name = cell.split("/")[-1].split(".")[0].split("_binsize_")[0]

		print(f"{cell_name}")
		cool_path = scool_path + "::" + cell
		cell_cool = Cooler(cool_path)
		cell_coverage = get_hic_coverage(cell_cool, binsize=COVERAGE_BINSIZE, excluded_chromosomes=EXCLUDED_CHROMS)
		cell_cnvs = call_cnvs(cell_coverage, binsize=COVERAGE_BINSIZE, ratio_col = "chromosomeNormalized_coverage", gaps=GAPS)
		cell_cnvs['is_significant'] = False
		cell_cnvs.loc[(cell_cnvs.genomewideNormalized_coverage - 1).abs() >= MAX_DEVIATION_FROM_GENOMEWIDE_MEDIAN, "is_significant"] = True
		cell_cnvs['cell'] = cell_name
		all_cell_cnvs.append(cell_cnvs)
	all_cell_cnvs = pd.concat(all_cell_cnvs, axis=0, ignore_index=True)

	def __get_cnv_type(x):
		if x.is_significant:
			if x['genomewideNormalized_coverage'] > 1:
				return "gain"
			elif x['genomewideNormalized_coverage'] < 1:
				return "loss"
			else:
				raise ValueError("Incompatible cnv type")
		else:
			return "NONE"

	all_cell_cnvs['cnv_type'] = all_cell_cnvs.apply(__get_cnv_type, axis=1)
	all_cell_cnvs = all_cell_cnvs[["cell", "chrom", "start", "end", "mean_value", "n_points", 
				   "coverage", "genomewideNormalized_coverage", "chromosomeNormalized_coverage", 
				   "is_significant", "cnv_type"]]
	all_cell_cnvs = all_cell_cnvs.rename(columns = {'chrom': 'chr'})
	all_cell_cnvs.to_csv(name, sep="\t", index=False, header=True)
else:
	all_cell_cnvs = pd.read_csv(name, sep="\t")


genome_size = Cooler(scool_path + "::" + list_coolers(scool_path)[0]).chromsizes.drop(EXCLUDED_CHROMS, axis=0).sum()
all_cells_cnv_covered_genome = all_cell_cnvs[all_cell_cnvs.is_significant]\
			 .assign(cnv_size = lambda x: x.end - x.start,
					 cnv_perc_genome = lambda x: x.cnv_size / genome_size,
					 n_cnvs = 1)\
			 .groupby(['cell'])[['cnv_perc_genome', 'n_cnvs']]\
			 .sum()\
			 .reindex(sorted(all_cell_cnvs.cell.unique()))\
			 .fillna(0)\
			 .reset_index()\
			 .assign(n_cnvs = lambda x: x.n_cnvs.astype(int))\
			 .merge(schic_samples.rename(columns = {'sample_name': 'cell'}), on = 'cell')\
			 .assign(condition = lambda x: x.cell_line + "_" + x.condition)\
			 .assign(cell_category = lambda x: x.apply(to_cell_category, axis=1))

all_cell_compartment_scores_with_cnvs = all_cell_compartment_scores[['cell', 'silhouette_score_comp']]\
											.merge(all_cells_cnv_covered_genome[['cell', 'cnv_perc_genome', 'n_cnvs', 'condition', 'cell_category', 'flag']], 
													on = 'cell')
selected_cells = all_cell_compartment_scores_with_cnvs[all_cell_compartment_scores_with_cnvs.flag.isna()]
palette = selected_cells[['cell', 'cell_category']]\
				.drop_duplicates()\
				.set_index('cell')\
				['cell_category']\
				.map(CELL_CATEGORY_PALETTE)\
				.to_dict()
ordered_cells = selected_cells[['cell', 'cnv_perc_genome']].drop_duplicates().sort_values("cnv_perc_genome")
cell_order = ordered_cells.cell.tolist()

bins = BedTool().window_maker(w = COVERAGE_BINSIZE, genome='hg19')
all_cell_cnvs_binned = bins\
			.intersect(BedTool.from_dataframe(all_cell_cnvs.loc[all_cell_cnvs.cell.isin(schic_samples[schic_samples.flag.isna()].sample_name),['chr', 'start', 'end', 'genomewideNormalized_coverage', 'cell']]).sort(), wa=True, wb=True)\
			.to_dataframe(names = ['chr', 'start', 'end', 'cnv_chr', 'cnv_start', 'cnv_end', 'genomewideNormalized_coverage', 'cell'])
bins = bins.to_dataframe(names = ['chr', 'start', 'end'])
bins = bins[bins.chr.isin(all_cell_cnvs.chr)]
chromosomes = [c for c in Cooler(scool_path + "::" + list_coolers(scool_path)[0]).chromnames if c not in EXCLUDED_CHROMS]
bins['chr'] = pd.Categorical(bins.chr, categories=chromosomes, ordered=True)
bins = bins.sort_values(['chr', 'start', 'end']).reset_index(drop=True)

limits = bins.assign(id = lambda x: x.index)\
			 .groupby('chr')['id'].agg(['min', 'max', "mean"])

heatmap = bins.merge(
				pd.pivot_table(all_cell_cnvs_binned, 
							   index = ['chr', 'start', 'end'], 
							   columns = "cell", 
							   values = "genomewideNormalized_coverage", 
							   fill_value=1)\
					.reset_index(), 
				on = ['chr', 'start', 'end'], how = "left")\
				.fillna(1)\
				.set_index(['chr', 'start', 'end'])\
				.reindex(index = bins)\
				.T\
				.reindex(index=cell_order)


g = plot_cnv_heatmap(heatmap, all_cells_cnv_covered_genome, limits)
g.savefig(os.path.join(OUTPUT_PATH, "all_CNVs_heatmap.pdf"))
plt.close(g.fig)


