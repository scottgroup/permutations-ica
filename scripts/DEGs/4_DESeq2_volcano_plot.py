import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v


# All colors
colors = {
    # Blue
    "testis": "#1f77b4ff", "cutadapt": "#1f77b4ff", "ensembl92": "#1f77b4ff",
    "HISAT2": "#1f77b4ff", "cufflinks": "#1f77b4ff",
    # Green
    "heart": "#2ca02cff", "refseq": "#2ca02cff",
    "tophat2": "#2ca02cff", "htseq": "#2ca02cff",
    # Orange
    "thyroid": "#ff7f0eff", "trimmomatic": "#ff7f0eff",
    "ensembl98": "#ff7f0eff", "STAR": "#ff7f0eff", "featureCounts": "#ff7f0eff",
    # Red
    "colon": "#da4144ff"
}


min_foldchange = snakemake.params.fold_change
min_pval = snakemake.params.pvalue


def get_pos_from_sample(sample):
    """ Calculate position in the grid """
    for i, (k, v) in enumerate(snakemake.config['tools'].items()):
        for tool in v:
            if tool in sample:
                return i


def plotting_scatter(ax, x, y, filt, color, cut=False):
    """ Creates the scatter plots """
    x = x[filt]
    y = y[filt]

    if cut:
        x = x[np.arange(0, len(x), 3)]
        y = y[np.arange(0, len(y), 3)]

    ax.scatter(x, y, color=color, alpha=1, s=2)


def plotting_nums(ax, ndown, nup):
    """ Adding number of genes on the volcano plots """
    ax.text(-25, 0, ndown)
    ax.text(25, 0, nup, horizontalalignment='right')


def plotting_volcano(ax, data, sample, min_pval, min_foldchange):
    """ Main filtering and plotting function """
    sample_down = sample.split('_')[0]
    sample_up = sample.split('_')[-1]

    # Data
    foldchange = data['log2FoldChange']
    pvalue = -np.log10(data['padj'])
    min_pval = -np.log10(min_pval)

    # Filters
    filt_up = (foldchange > min_foldchange) & (pvalue > min_pval)
    filt_down = (foldchange < -min_foldchange) & (pvalue > min_pval)
    filt_mass = (foldchange < 2) & (pvalue < 35)
    filt_other = ~(filt_up | filt_down | filt_mass)

    # Plotting
    plotting_scatter(ax, foldchange, pvalue, filt_up, colors[sample_up])
    plotting_scatter(ax, foldchange, pvalue, filt_down, colors[sample_down])
    plotting_scatter(ax, foldchange, pvalue, filt_mass, 'k', cut=True)
    plotting_scatter(ax, foldchange, pvalue, filt_other, 'k')

    plotting_nums(ax, np.sum(filt_down), np.sum(filt_up))

    ax.set_xlim(-32, 32)
    ax.set_ylim(-10, 320)
    ax.set_title(sample)


nrows = 4
ncols = 3
used_subplots = np.zeros((nrows, ncols))

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=mm2inch((84,150)))
x = dict()
for csv in snakemake.input.DEGs:
    sample = csv.split('.')[0].split('/')[-1]

    # Getting position in subplots
    y = get_pos_from_sample(sample)
    if y not in x:
        x[y] = 0
    else:
        x[y] += 1

    # Loading data
    data = pd.read_csv(csv, sep=',')
    # Correcting for p-value == 0, replacing with approx smallest value
    data.loc[data['padj'] == 0, 'padj'] = 2.5e-308

    plotting_volcano(axes[y, x[y]], data, sample, min_pval=min_pval, min_foldchange=min_foldchange)
    if x[y] != 0:
        axes[y, x[y]].set_yticks([])
    else:
        axes[y, x[y]].set_ylabel('-log(p-value)')
        axes[y, x[y]].set_xlabel('log2(FoldChange)')
    used_subplots[y, x[y]] = 1


# Formatting graph
for x in range(used_subplots.shape[0]):
    for y in range(used_subplots.shape[1]):
        if used_subplots[x, y] == 0:
            axes[x, y].axis('off')

plt.subplots_adjust(wspace=0, hspace=0.45, right=0.95)
plt.savefig(snakemake.output.plot, dpi=1000)
