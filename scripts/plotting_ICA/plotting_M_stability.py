import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v


# Loading data
data = pd.read_csv(
    snakemake.input.data, names=['IDs', 'metric'], sep='\t', header=0
)

# Half IDs
IDs = [id[1:] if idx%2 == 0 else '' for idx, id in enumerate(data['IDs'].tolist())]

# Params
max = int(snakemake.wildcards.max)
min = int(snakemake.wildcards.min)

plt.figure(figsize=mm2inch(78,50))
plt.plot(data['IDs'], data['metric'], 'k', linewidth=2)
plt.xlim(min, max-min)
locs, labels = plt.xticks()
plt.xticks(locs, IDs)

plt.ylabel('Mean Squared Error')
plt.xlabel('Number of expression modes')
plt.tight_layout()
plt.savefig(snakemake.output.plot)
