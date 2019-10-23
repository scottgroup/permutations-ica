import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



def draw_scatter(ax, df, col1, col2, cols, var='tissue'):

    # Get the different values for a variable
    var_list = df.index.get_level_values(var)
    vars = list(set(var_list))

    for var in vars:
        filt = var_list == var
        x = df[cols[col1]][filt]
        y = df[cols[col2]][filt]

        ax.scatter(
            x, y,
            s=70,
            alpha=0.35, edgecolor='k',
            label=var
        )

# Importing projection
proj = pd.read_csv(snakemake.input.proj, sep='\t', index_col=[0,1,2,3,4,5])

# Importing components variable
comp_var = pd.read_csv(snakemake.input.comp_list, sep='\t', names=['comp', 'variable'])

# Heart : 0, 2, 7, 12
# Colon : 4, 6, 8, 13
# Testis : 3
# Thyroid : 10, 11, 14, 15
# TODO : automatic selection of components
cols = ['0', '3', '13', '15']
n_grid = len(cols)

fig, axes = plt.subplots(n_grid, n_grid)

for i in range(n_grid):
    for j in range(i+1):
        draw_scatter(axes[i, j], proj, i, j, cols)
        # axes[i,j].set_xlim(axes[i,j].get_xlim()[1], axes[i,j].get_xlim()[0])

    axes[n_grid-1,j].set_xlabel('Comp ' + cols[i])
    axes[i,0].set_ylabel('Comp ' + cols[i])

# Formatting
for i in range(n_grid):
    for j in range(n_grid):
        axes[i,j].set_xticks([])
        axes[i,j].set_yticks([])

        if j > i:
            axes[i, j].axis('off')


handles, labels = axes[0,0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.85,0.85))

plt.subplots_adjust(hspace=0, wspace=0)
plt.show()
