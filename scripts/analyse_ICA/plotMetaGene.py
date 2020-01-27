import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v


pos = {
    '5': {
        'xticks': [0, 10, 90],
        'xticklabels': ['-10', '0', '+80']
    },
    '3': {
        'xticks': [90, 80, 0],
        'xticklabels': ['-10', '0', '+80']
    }
}

# Loading data
mean = pd.read_csv(snakemake.input.mean, sep='\t', header=[0,1], index_col=0)
std = pd.read_csv(snakemake.input.std, sep='\t', header=[0,1])

# Plotting
figure, axes = plt.subplots(nrows=1, ncols=2, figsize=mm2inch(178/7*2.5, 50))

x = np.arange(len(mean))
handles = dict()
for ax, prime in enumerate(['5', '3']):

    for tool in mean.columns.get_level_values('tool'):
        quant = mean.loc(axis=1)[prime, tool]
        std_quant = std.loc(axis=1)[prime, tool]

        handles[tool], = axes[ax].plot(
            x, quant,
            color=colors[tool], linewidth=2
        )
        handles[tool].set_label(tool)

    # Setting limits
    xmin, xmax = 0, len(x)-1
    ymin, ymax = 0, 1
    axes[ax].set_xlim([xmin, xmax])
    axes[ax].set_ylim([ymin, ymax])

    # Removing axis and plotting xline
    axes[ax].set_frame_on(False)
    axes[ax].get_xaxis().tick_bottom()
    axes[ax].axes.get_yaxis().set_visible(False)
    axes[ax].add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=1))

    # Drawing vertical lines
    axes[ax].axvline(pos[prime]['xticks'][1], linestyle=':', color='k')

    # Setting axes labels
    axes[ax].set_xticks(pos[prime]['xticks'])
    axes[ax].set_xticklabels(pos[prime]['xticklabels'])

axes[0].legend(handles.values(), labels=['HISAT2', 'STAR', 'TopHat2'])

figure.tight_layout()
plt.savefig(snakemake.output.plot)
