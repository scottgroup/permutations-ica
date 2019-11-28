import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

k_neigh = int(snakemake.params.k_neighbour)
vars = ['tissue', 'trimmer', 'annotation', 'aligner', 'quantifier']

# Importing projection
proj = pd.read_csv(snakemake.input.proj, sep='\t', index_col=[0,1,2,3,4,5])


# Initialize score dataframe
index = list()
tool_dict = dict()
for k, v in zip(proj.index.names, proj.index.levels):
    tool_dict[k] = list(v)
for var in vars:
    index.extend(tool_dict[var])
score_df = pd.DataFrame(index=index, columns=proj.columns, dtype=float)

# Need to compute k-neighbour identity
for c in proj.columns:
    for var in vars:
        col_series = proj[c].copy(deep=True)
        col_series.index = col_series.index.get_level_values(var)

        # Col-wise difference matrix
        repmat = np.tile(col_series.values, (len(col_series), 1)).T
        diff_mat = np.abs(repmat - col_series.values)
        score = list()
        for row in range(diff_mat.shape[0]):
            idxs = np.argpartition(diff_mat[row, :], k_neigh)[:k_neigh]
            identity = col_series.index[row]
            score.append(
                np.sum(col_series.iloc[idxs].index == identity) / k_neigh
            )
        score_series = pd.Series(score, index=col_series.index)
        resp_dict = score_series.groupby(var).mean().to_dict()

        for k, v in resp_dict.items():
            score_df.loc[k, c] = v


# Plotting
height = list()
for var in vars:
    height.append(len(tool_dict[var]))

fig, axes = plt.subplots(
    nrows=len(vars), ncols=1,
    gridspec_kw={'height_ratios': height},
    figsize=(14, 12)
)

for i, var in enumerate(vars):
    sns.heatmap(
        score_df[score_df.index.isin(tool_dict[var])],
        linewidths=0.5,
        annot=True, ax=axes[i],
        cmap="YlGnBu", vmax=1, vmin=0.30
    )
    axes[i].set_yticklabels(labels=axes[i].get_yticklabels(), rotation=0)

    if i < len(vars)-1:
        axes[i].set_xticks([])
    else:
        axes[i].set_xlabel('Components')


plt.subplots_adjust(hspace=0.1, top=0.95, right=0.95, bottom=0.1)
plt.savefig(snakemake.output.plot)
