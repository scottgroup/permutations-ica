import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import permutations

from scipy.cluster.hierarchy import dendrogram, linkage, distance, fcluster
from scipy.spatial import distance_matrix
from sklearn.cluster import AgglomerativeClustering


def get_bdiagonal(mat, threshold=0.9):
    """
    """
    d = np.shape(mat)[0]
    blocks = np.zeros((d, d))

    for i in range(int(d/n)):
        low = i*n
        high = i*n + n
        blocks[low:high, low:high] = 1

    return blocks

def get_error(mat, bdiagonal):
    """ Calculating the MSE """
    err = np.mean((mat - bdiagonal)**2)
    return err


def getting_components(path):
    # Reading data
    components = pd.read_csv(
        path,
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
    return X


# Plotting
IDs = list()
metrics = list()

data = list()
n = int(snakemake.wildcards.n)
for file in snakemake.input:
    mat = getting_components(file)
    bdiagonal = get_bdiagonal(mat)
    metric = get_error(mat, bdiagonal)
    M_id = 'M' + file.split('/M')[-1].split('_')[0]
    data.append([M_id, metric])

df = pd.DataFrame(data)
print(df)
df.to_csv(snakemake.output.data, sep='\t', index=None)
