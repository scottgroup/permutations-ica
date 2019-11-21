import numpy as np


def logcosh_gs(s, a=1):
    """
    Calculating the negentropy using the logcosh approximation
    """
    return np.mean(
        np.log(np.cosh(a * s))/a
    )


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
    
    Code from :
    https://stackoverflow.com/questions/31528800/how-to-implement-zca-whitening-python

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
