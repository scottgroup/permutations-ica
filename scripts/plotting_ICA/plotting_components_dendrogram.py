import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy

import sys
sys.setrecursionlimit(10000)

# Loading data
data = pd.read_csv(snakemake.input.components, sep='\t', index_col=[0,1,2,3])
threshold = snakemake.params.threshold

# Calculating hierarchy
Z = hierarchy.linkage(data.T, 'single')

# Graphing
plt.figure(figsize=(20,20))
dn = hierarchy.dendrogram(
    Z,
    labels=data.columns,
    color_threshold=(1-threshold)*max(Z[:,2])
)

# Saving to file
plt.savefig(snakemake.output.plot)
plt.savefig(snakemake.output.plot2)
