import numpy as np
import pandas as pd

from Dataset import Dataset
from pyProDenICA import pyProDenICA


def run_ICA(dataset, M, max_it):
    """
    Run the ICA from a model and outputs components in a pandas
    DataFrame
    """

    data = np.array(dataset.data).T

    std = np.std(data, axis=1)
    mean_std = np.mean(std)
    data = data[std > mean_std*2, :]

    print('Shape data')
    print(data.shape)


    W0, crit0, nw, n_it = pyProDenICA(
        X=data,
        k=M, maxit=max_it,
        whiten=True,
        trace=True
    )

    components_i = pd.DataFrame(
        W0,
        columns=dataset.data.columns
    )
    components_i = components_i.transpose()
    components_i.columns = [
        'r_' + str(i) + ' n_' + str(c) for c in components_i.columns
    ]

    return components_i, crit0, n_it


# Importating, formating and scaling the data
data_slice = snakemake.config['ICA_datasets'][snakemake.wildcards.dataset]
dataset = Dataset(snakemake.input.dataset, data_slice)

# Setting parameters
M = int(snakemake.wildcards.M)  # Number of components
n = int(snakemake.wildcards.n)  # Number of bootstrapping iteration
max_it = int(snakemake.params.max_it)

# Bootstrapping the ICA
fit_min = list()
for i in range(n):
    _components, crit0, n_it = run_ICA(dataset, M, max_it)

    if n_it < max_it:
        fit_min.append(crit0)
        if i == 0:
            components = _components
        else:
            components = pd.concat([components, _components], axis=1)

# Writing components to disk
components.to_csv(snakemake.output.components, sep='\t')

# Writing fit minimum to disk
with open(snakemake.output.fit_min, 'w') as f:
    for min in fit_min:
        f.write(min + '\n')
