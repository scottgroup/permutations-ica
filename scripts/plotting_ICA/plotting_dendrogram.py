import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy

import sys
sys.setrecursionlimit(10000)

# Loading data
data = pd.read_csv(snakemake.input.components, sep='\t', index_col=[0,1,2,3])
print(data)

Z = hierarchy.linkage(data.T, 'single')
plt.figure()
dn = hierarchy.dendrogram(Z)

plt.savefig(snakemake.output.plot)
plt.savefig(snakemake.output.plot2)
