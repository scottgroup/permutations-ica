import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

min_foldchange = snakemake.params.fold_change
min_pval = snakemake.params.pvalue


def get_pos_from_sample(sample):
    """ """
    tissues = snakemake.config['datasets'].keys()
    for tissue in tissues:
        if tissue in sample:
            return 0

    for i, (k, v) in enumerate(snakemake.config['tools'].items()):
        for tool in v:
            if tool in sample:
                return i+1


def plotting_scatter(ax, x, y, filt, color):
    """ """
    ax.scatter(
        x[filt], y[filt],
        color=color, alpha=0.2
    )


def plotting_nums(ax, ndown, nup):
    """ """
    ax.text(-25, 0, ndown)
    ax.text(25, 0, nup, horizontalalignment='right')


def plotting_volcano(ax, data, sample, min_pval=0.0001, min_foldchange=2):
    """ """
    # Data
    foldchange = data['log2FoldChange']
    pvalue = -np.log10(data['padj'])

    # Filters
    filt_up = (foldchange > min_foldchange) & (pvalue > min_pval)
    filt_down = (foldchange < -min_foldchange) & (pvalue > min_pval)
    filt_other = ~(filt_up | filt_down)

    # Plotting
    plotting_scatter(ax, foldchange, pvalue, filt_up, 'g')
    plotting_scatter(ax, foldchange, pvalue, filt_down, 'r')
    plotting_scatter(ax, foldchange, pvalue, filt_other, 'k')

    plotting_nums(ax, np.sum(filt_down), np.sum(filt_up))

    ax.set_xlim(-31, 31)
    ax.set_ylim(-25, 325)
    ax.set_title(sample)




nrows = 5
ncols = 6
used_subplots = np.zeros((nrows, ncols))

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20,12))
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
    if (y != 0) | (x[y] != 0):
        axes[y, x[y]].set_xticks([])
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

plt.suptitle(snakemake.wildcards.dataset)
plt.subplots_adjust(wspace=0)

for x in range(used_subplots.shape[0]):
    for y in range(used_subplots.shape[1]):
        if used_subplots[x, y] == 1:
            num_missing = ncols - np.sum(used_subplots[x, :])
            box = axes[x,y].get_position()
            offset = (box.x1 - box.x0) * num_missing / 2
            box.x0 = box.x0 + offset
            box.x1 = box.x1 + offset
            axes[x,y].set_position(box)

plt.savefig(snakemake.output.plot)
