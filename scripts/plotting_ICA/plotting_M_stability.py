import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import permutations

from scipy.cluster.hierarchy import dendrogram, linkage, distance, fcluster
from scipy.spatial import distance_matrix
from sklearn.cluster import AgglomerativeClustering


def finding_blocks(mat, threshold=0.9):
    """
    """
    d = np.shape(mat)[0]
    blocks = np.zeros(d)

    # Walking the blocks
    i = 0
    step = 0
    block_num = 0
    while i + step < d:
        sub_mat = mat[i:i+1+step, i:i+1+step]
        blocks[i + step] = block_num
        if np.min(sub_mat) < threshold:
            i += step

            # Incrementing if diagonal is lower then threshold
            if step == 0:
                i += 1

            step = 0
            block_num += 1
        else:
            step += 1

    return blocks


def pearson_correlation(mat, blocks=[]):
    """
    Adapted the idea from :
    https://math.stackexchange.com/questions/1392491/measure-of-how-much-diagonal-a-matrix-is

    """
    # Setting var
    d = np.shape(mat)[0]
    j = np.ones(d)
    r = np.arange(1, d+1, 1)

    xmat = np.outer(r, j.T)
    # If using the block scheme
    if len(blocks) > 1:
        for i in range(int(np.max(blocks)+1)):
            args = np.argwhere(blocks==i)
            for x,y in permutations(args, 2):
                xmat[x, y] = (xmat[x, x] + xmat[y, y])/2
    ymat = xmat.T

    n = d*d
    sum_x = np.sum( xmat * mat )
    sum_x2 = np.sum( xmat**2 * mat )
    sum_y = np.sum( ymat * mat )
    sum_y2 = np.sum( ymat**2 * mat )
    sum_xy = np.sum( xmat * ymat * mat )

    num = n * sum_xy - sum_x * sum_y
    dem = np.sqrt(n*sum_x2 - sum_x**2) * np.sqrt(n*sum_y2 - sum_y**2)

    return num/dem


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
pearsons = list()

for file in snakemake.input:
    mat = getting_components(file)
    blocks = finding_blocks(mat)
    pearson = pearson_correlation(np.abs(mat), blocks)
    M_id = 'M' + file.split('/M')[-1].split('_')[0]
    IDs.append(M_id)
    pearsons.append(pearson)

plt.plot(IDs, pearsons)
plt.savefig(snakemake.output.plot)
