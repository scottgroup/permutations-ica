import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v


def draw_scatter(ax, df, col1, col2, cols, var='tissue'):
    """ Main plotting fun """
    # Get the different values for a variable
    var_list = df.index.get_level_values(var)
    vars = list(set(var_list))

    for var in vars:
        filt = var_list == var
        x = df[cols[col1]][filt]
        y = df[cols[col2]][filt]

        ax.scatter(
            x, y,
            s=25,
            alpha=0.35, edgecolor='k',
            label=var
        )

# Importing projection
proj = pd.read_csv(snakemake.input.proj, sep='\t', index_col=[0,1,2,3,4,5])

# Importing components variable
comp_var = pd.read_csv(snakemake.input.comp_list, sep='\t', names=['comp', 'variable'])
comps = snakemake.params.comps
n_grid = len(comps)

fig, axes = plt.subplots(n_grid, n_grid, figsize=mm2inch((86,86)))

for i in range(n_grid):
    for j in range(i+1):
        draw_scatter(axes[i, j], proj, i, j, comps)

    axes[n_grid-1,j].set_xlabel('Mode ' + str(int(comps[i]) + 1))
    axes[i,0].set_ylabel('Mode ' + str(int(comps[i]) + 1))

# Formatting
for i in range(n_grid):
    for j in range(n_grid):
        axes[i,j].set_xticks([])
        axes[i,j].set_yticks([])

        # Modifying lims
        scale = 0.025
        left, right = axes[i,j].get_xlim()
        axes[i,j].set_xlim((left-scale, right+scale))
        bottom, top = axes[i,j].get_ylim()
        axes[i,j].set_ylim((bottom-scale, top+scale))

        if j > i:
            axes[i, j].axis('off')


handles, labels = axes[0,0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.85,0.85))
plt.subplots_adjust(hspace=0, wspace=0)
plt.savefig(snakemake.output.plot)
