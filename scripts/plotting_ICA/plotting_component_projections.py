import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

# File path
fpath = snakemake.params.fpath.format(**snakemake.wildcards)
os.makedirs(os.path.dirname(fpath), exist_ok=True)

# Different entries in the plot
categories = ['tissue', 'trimmer', 'annotation', 'aligner', 'quantifier']

# Setting up the figure
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 12

# Input projection data
projection = pd.read_csv(
    snakemake.input.projection, sep='\t', index_col=[0, 1, 2, 3, 4, 5]
)

# Plotting for each variable
for comp in projection.columns.tolist():
    fig, axes = plt.subplots(len(categories), figsize=(16, 12))

    # Creating Y jitter position
    y = np.random.rand(len(projection))

    for category, idx in zip(categories, range(len(categories))):
        data_dict = projection[comp].groupby(category).apply(list).to_dict()

        for key, val in data_dict.items():
            filt = projection[comp].index.get_level_values(category) == key
            slice = np.arange(len(filt))

            axes[idx].scatter(
                val, y[slice[filt]],
                label=key,
                alpha=0.4,
                s=75
            )
            axes[idx].legend(
                bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                ncol=len(key), mode="expand", borderaxespad=0.
           )
            axes[idx].set_xticks([], [])
            axes[idx].set_yticks([], [])
            axes[idx].set_ylim([-0.2, 1.2])

    plt.subplots_adjust(wspace=0, hspace=0.325)
    plt.savefig(fpath.format(comp=comp))
