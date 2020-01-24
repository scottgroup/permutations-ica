import numpy as np
import pandas as pd

from Dataset import Dataset
from sklearnFastICA import sklearnFastICA
import running_ICA_utils as rICA


def run_ICA(data, M, max_it, tolerance, model_i):
    """
    Run the ICA from a model and outputs components in a pandas
    DataFrame
    """

    components = sklearnFastICA(data, M, max_it, tolerance)

    # If convergence
    if components is not None:
        # Calculting the neg-entropy
        obj_f = rICA.logcosh_gs(components)

        # Formatting the components
        components_i = pd.DataFrame(components, index=data.index)
        components_i.columns = [
            'r_' + str(model_i) + ' n_' + str(c) for c in components_i.columns
        ]
        return components_i, obj_f

    # If not converged
    else:
        return None, None


# Importating, formating and scaling the data
data_slice = snakemake.config['ICA_datasets'][snakemake.wildcards.dataset]['variables']
dataset = Dataset(snakemake.input.counts, data_slice)

# Setting parameters
M = int(snakemake.wildcards.M)  # Number of components
n = int(snakemake.wildcards.n)  # Number of bootstrapping iteration
max_it = int(snakemake.params.max_it)
tolerance = float(snakemake.params.tolerance)
std_from_mean = float(snakemake.wildcards.std)

# Preparing data
data = rICA.prepare_data(dataset, std_from_mean)

# Whitening the data
data = pd.DataFrame(
    rICA.zca_whitening_matrix(data),
    index=data.index, columns=data.columns
)

# Iterating the ICA
fit_min = list()
i = 0
while i < n:
    # Subselecting data if bootstrapped
    if 'boot' in snakemake.wildcards.keys():
        n_columns = len(data.columns)
        rand_columns = np.random.choice(
            n_columns,
            int(int(snakemake.wildcards.boot)/100*n_columns),
            replace=False
        )
        sub_data = data[data.columns[rand_columns]]
    else:
        sub_data = data

    _components, crit0 = run_ICA(sub_data, M, max_it, tolerance, i)

    if _components is not None:
        fit_min.append(crit0)
        if i == 0:
            components = _components
        else:
            components = pd.concat([components, _components], axis=1)
        i += 1


# Writing components to disk
components.to_csv(snakemake.output.components, sep='\t')

# Writing fit minimum to disk
with open(snakemake.output.fit_min, 'w') as f:
    for min in fit_min:
        f.write(str(min) + '\n')
