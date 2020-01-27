import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v


def reading_file(path):
    up, down = list(), list()
    with open(path, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line == '>Positive genes':
                reading_up = True
            elif line == '>Negative genes':
                reading_up = False
            else:
                if reading_up:
                    up.append(line)
                else:
                    down.append(line)
    return up, down


def get_pos_from_sample(sample):
    """ Calculate position in the grid """
    for i, (k, v) in enumerate(snakemake.config['tools'].items()):
        for tool in v:
            if tool in sample:
                return i


def plotting_scatter(ax, x, y, color='k'):
    """ Creates the scatter plots """
    ax.scatter(x, y, alpha=1, s=2, color=color)


def plotting_nums(ax, ndown, nup):
    """ Adding number of genes on the volcano plots """
    ax.text(-25, 0, ndown)
    ax.text(25, 0, nup, horizontalalignment='right')


def plotting_volcano(ax, data, sample):
    """ Main filtering and plotting function """
    # Data
    data['padj'] = -np.log10(data['padj'])

    # Data up and down
    up_data = data[data[data.columns[0]].isin(
        [ens2symbol[g] for g in up]
    )]
    down_data = data[data[data.columns[0]].isin(
        [ens2symbol[g] for g in down]
    )]


    plotting_scatter(ax, data['log2FoldChange'], data['padj'])
    plotting_scatter(ax, up_data['log2FoldChange'], up_data['padj'], 'r')
    plotting_scatter(ax, down_data['log2FoldChange'], down_data['padj'], 'b')

    ax.set_xlim(-32, 32)
    ax.set_ylim(-10, 320)
    ax.set_title(sample)


# Reading component file
fcomp = snakemake.input.gene_list
up, down = reading_file(fcomp)

# Loading HGNC data
hgnc_df = pd.read_csv(snakemake.input.HGNC, sep='\t')
ens2symbol = hgnc_df.set_index('ensembl_id')['symbol'].to_dict()

# Getting into drawing
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

    plotting_volcano(axes[y, x[y]], data, sample)
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
