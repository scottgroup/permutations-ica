import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v

# Names
names = {
    "colon": "Colon", "heart": "Heart", "testis": "Testis", "thyroid": "Thyroid",
    "cutadapt": "Cutadapt", "trimmomatic": "Trimmomatic",
    "ensembl92": "Ensembl92", "ensembl98": "Ensembl98", "refseq": "RefSeq",
    "tophat2": "TopHat2", "HISAT2": "HISAT2", "STAR": "STAR",
    "cufflinks": "Cufflinks", "featureCounts": "featureCounts", "htseq": "HTSeq"
}

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

# File path
fpath = snakemake.params.fpath.format(**snakemake.wildcards)
os.makedirs(os.path.dirname(fpath), exist_ok=True)

# Different entries in the plot
categories = ['tissue', 'trimmer', 'annotation', 'aligner', 'quantifier']

# Input projection data
projection = pd.read_csv(
    snakemake.input.projection, sep='\t', index_col=[0, 1, 2, 3, 4, 5]
)

# Plotting for each variable
for comp in projection.columns.tolist():
    fig, axes = plt.subplots(len(categories), figsize=mm2inch((178, 115)), sharex=True)
    # Creating Y jitter position
    y = np.random.rand(len(projection))

    for category, idx in zip(categories, range(len(categories))):
        data_dict = projection[comp].groupby(category).apply(list).to_dict()
        for key, val in data_dict.items():
            filt = projection[comp].index.get_level_values(category) == key
            slice = np.arange(len(filt))

            sns.distplot(
                val, rug=False, hist=False,
                kde=True,
                kde_kws={"lw": 3, "label": names[key], "cut": 2, "bw":0.01},
                rug_kws={"height": 0.15, "alpha":0.3},
                color=colors[key],
                ax=axes[idx]
            )

        axes[idx].legend(loc="center right")
        axes[idx].set_yticks([], [])
        axes[idx].set_ylabel(category.capitalize())
        axes[idx].yaxis.labelpad = 0
    axes[idx].set_xlabel('Projection')

    xlim = axes[idx].get_xlim()
    axes[idx].set_xlim((xlim[0] + (xlim[1]-xlim[0])*0.025, xlim[1] + (xlim[1]-xlim[0])*0.075))

    plt.suptitle('Expression mode ' + str(int(comp)+1))
    plt.subplots_adjust(wspace=0, hspace=0, top=0.95, bottom=0.05, left=0.025, right=0.975)
    plt.savefig(fpath.format(comp=comp))
