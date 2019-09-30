import numpy as np
import pandas as pd

from Dataset import Dataset
from ICAmethods.sklearnFastICA import sklearnFastICA
# from pyProDenICA import pyProDenICA, ProDenICA
# from pyFastICA import running_fastICA


def prepare_data(dataset, std_from_mean):
    """
    Centering data, and keeping only genes that are at least std_from_mean std
    from the std mean.

    Use std_from_mean = 0 to keep all genes.
    """
    data_std = dataset.data.std(axis=1)
    mean_std = np.mean(data_std)
    data = dataset.data[data_std >= std_from_mean*mean_std]

    return data


def zca_whitening_matrix(X):
    """
    Function to compute ZCA whitening matrix (aka Mahalanobis whitening).
    INPUT:  X: [M x N] matrix.
        Rows: Variables
        Columns: Observations
    OUTPUT: ZCAMatrix: [M x M] matrix -> Not true anymore
    """
    # Covariance matrix [column-wise variables]: Sigma = (X-mu)' * (X-mu) / N
    sigma = np.cov(X, rowvar=True)  # [M x M]
    # Singular Value Decomposition. X = U * np.diag(S) * V
    U, S, V = np.linalg.svd(sigma)
        # U: [M x M] eigenvectors of sigma.
        # S: [M x 1] eigenvalues of sigma.
        # V: [M x M] transpose of U
    # Whitening constant: prevents division by zero
    epsilon = 1e-5
    # ZCA Whitening matrix: U * Lambda * U'
    zca_whitener = np.dot(U, np.dot(np.diag(1.0/np.sqrt(S + epsilon)), U.T))  # [M x M]
    return np.dot(zca_whitener, X)


def logcosh_gs(s, a=1):
    """
    Calculating the negentropy using the logcosh approximation
    """
    return np.mean(
        np.log(np.cosh(a * s))/a
    )


def run_ICA(ICAmethod, data, M, max_it, tolerance, model_i):
    """
    Run the ICA from a model and outputs components in a pandas
    DataFrame
    """

    if ICAmethod == 'rProDenICA':
        # Using rProDenICA
        pass

    elif ICAmethod == 'pyProDenICA':
        # Using pyProDenICA
        pass

    elif ICAmethod == 'sklearnFastICA':
        # Using sklearnFastICA
        components = sklearnFastICA(data, M, max_it, tolerance)

    elif ICA_method == 'cuFastICA':
        # Using customFastICA
        pass

    # Using fastICA
    # W, obj_f, i = running_fastICA(np.array(data), M)
    # components = W

    # Using pyProDenICA
    # W0, components, obj_f, nw, n_it = ProDenICA(
    #     X=np.array(data),
    #     k=M, maxit=max_it,
    #     whiten=True,
    #     trace=True,
    #     restarts=10
    # )

    # If convergence
    if components is not None:
        # Calculting the neg-entropy
        obj_f = logcosh_gs(components)

        # Formatting the components

        # print(components[:13, :])

        components_i = pd.DataFrame(components, index=data.index)

        # print(components_i)

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
ICAmethod = snakemake.wildcards.ICAmethod
M = int(snakemake.wildcards.M)  # Number of components
n = int(snakemake.wildcards.n)  # Number of bootstrapping iteration
max_it = int(snakemake.params.max_it)
tolerance = float(snakemake.params.tolerance)
std_from_mean = float(snakemake.wildcards.std)

# Preparing data
data = prepare_data(dataset, std_from_mean)

# Whitening the data
data = pd.DataFrame(
    zca_whitening_matrix(data),
    index=data.index, columns=data.columns
)

# Bootstrapping the ICA
fit_min = list()
i = 0
while i < n:
    _components, crit0 = run_ICA(ICAmethod, data, M, max_it, tolerance, i)

    if _components is not None:
        print('Success ', str(i))
        fit_min.append(crit0)
        if i == 0:
            components = _components
        else:
            components = pd.concat([components, _components], axis=1)
        i += 1
    else:
        print('\t Retrying ', str(i))


# Writing components to disk
components.to_csv(snakemake.output.raw_components, sep='\t')

# Writing fit minimum to disk
with open(snakemake.output.fit_min, 'w') as f:
    for min in fit_min:
        f.write(str(min) + '\n')
