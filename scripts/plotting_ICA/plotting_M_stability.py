import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')

# Loading data
data = pd.read_csv(
    snakemake.input.data, names=['IDs', 'metric'],
    sep='\t', header=0
)

# Half IDs
IDs = [id[1:] if idx%2 == 0 else '' for idx, id in enumerate(data['IDs'].tolist())]
print(IDs)

plt.plot(data['IDs'], data['metric'], 'k', linewidth=2)
plt.xlim(6, 29)
locs, labels = plt.xticks()
plt.xticks(locs, IDs)

plt.ylabel('Mean Squared Error')
plt.xlabel('Number of expression modes')
plt.savefig(snakemake.output.plot)
