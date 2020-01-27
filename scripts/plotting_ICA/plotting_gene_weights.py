import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import snakemake as smk

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v

# Read data
data = pd.read_csv(
    snakemake.input.components_mean, sep='\t', index_col=[0,1,2,3]
)

max = np.max(np.max(data))*1.1
min = np.min(np.min(data))*1.1

for c in data.columns.tolist():
    weights = sorted(data[c].values)
    x = np.arange(len(weights))

    plt.figure(figsize=(7,6))
    plt.plot(x, weights, 'k')
    plt.ylim((min,max))
    plt.xlabel('Genes')
    plt.ylabel('Weight')
    plt.title('Expression mode ' + str(c))

    fname = snakemake.params.dir + '/comp{c}.svg'.format(c=c)
    plt.savefig(fname)

smk.shell("touch {snakemake.output.tkn}")
