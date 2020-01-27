import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, distance, fcluster
from scipy.spatial import distance_matrix

from sklearn.cluster import AgglomerativeClustering

# Reading data
components = pd.read_csv(
    snakemake.input.components,
    sep='\t', index_col=[0, 1, 2, 3]
)
X = components.corr().values

"""
Inspired from
https://github.com/TheLoneNut/CorrelationMatrixClustering/blob/master/CorrelationMatrixClustering.ipynb
"""

cluster_th = 4

X = components.corr().values
d = distance.pdist(X)
L = linkage(d, method='complete')
ind = fcluster(L, 0.25*d.max(), 'distance')

columns = [components.columns.tolist()[i] for i in list(np.argsort(ind))]
components = components.reindex(columns, axis=1)

unique, counts = np.unique(ind, return_counts=True)
counts = dict(zip(unique, counts))

i = 0
j = 0
columns = []
for cluster_l1 in set(sorted(ind)):
    j += counts[cluster_l1]
    sub = components[components.columns.values[i:j]]
    if counts[cluster_l1]>cluster_th:
        X = sub.corr().values
        d = distance.pdist(X)
        L = linkage(d, method='complete')
        ind = fcluster(L, 0.25*d.max(), 'distance')
        col = [sub.columns.tolist()[i] for i in list((np.argsort(ind)))]
        sub = sub.reindex(col, axis=1)
    cols = sub.columns.tolist()
    columns.extend(cols)
    i = j
components = components.reindex(columns, axis=1)

# Get clusters
X = components.corr().values
d = distance.pdist(X)
L = linkage(d, method='complete')
ind = fcluster(L, 0.1*d.max(), 'distance')

# Creating the figure
fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111)
cax = ax.imshow(
    components.corr().values,
    interpolation=None, vmin=-1, vmax=1,
    cmap=mpl.cm.get_cmap('seismic_r')
)
# Adding colorbar
fig.colorbar(cax, norm=mpl.colors.Normalize(vmin=-1, vmax=1))

# Saving graph to file
plt.savefig(snakemake.output.plot)
