"""
https://github.com/akcarsten/Independent_Component_Analysis

"""
import numpy as np

def center(x):
    mean = np.mean(x, axis=1, keepdims=True)
    centered =  x - mean
    return centered, mean


def covariance(x):
    mean = np.mean(x, axis=1, keepdims=True)
    n = np.shape(x)[1] - 1
    m = x - mean

    return (m.dot(m.T))/n


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


def whiten(x):
    # Calculate the covariance matrix
    coVarM = covariance(x)

    # Single value decoposition
    U, S, V = np.linalg.svd(coVarM)

    # Calculate diagonal matrix of eigenvalues
    d = np.diag(1.0 / np.sqrt(S))

    # Calculate whitening matrix
    whiteM = np.dot(U, np.dot(d, U.T))

    # Project onto whitening matrix
    Xw = np.dot(whiteM, x)

    return Xw, whiteM


def fastICA(signals, M, alpha = 1, thresh=1e-11, iterations=5000):
    m, n = signals.shape
    k = M

    # Initialize random weights
    W = np.random.rand(m, k).T

    for c in range(k):
        w = W[c, :].copy().reshape(m, 1)
        w = w / np.sqrt((w ** 2).sum())

        i = 0
        lim = 100
        while ((lim > thresh) & (i < iterations)):

            # Dot product of weight and signal
            ws = np.dot(w.T, signals)

            # Pass w*s into contrast function g
            wg = np.tanh(ws * alpha).T

            # Pass w*s into g prime
            wg_ = (1 - np.square(np.tanh(ws))) * alpha

            # Update weights
            wNew = (signals * wg.T).mean(axis=1) - wg_.mean() * w.squeeze()

            # Decorrelate weights
            wNew = wNew - np.dot(np.dot(wNew, W[:c].T), W[:c])
            wNew = wNew / np.sqrt((wNew ** 2).sum())

            # Calculate limit condition
            lim = np.abs(np.abs((wNew * w).sum()) - 1)
            # print('\tIteration ', i, ' limit ', lim)

            # Update weights
            w = wNew

            # Update counter
            i += 1

        W[c, :] = w.T
    neg = np.mean(
        np.log(np.cosh(alpha*np.dot(W, signals)))/alpha
    )
    return W.T, neg, i


def running_fastICA(x, M):
    """ """
    # Center signals
    Xc, meanX = center(x)

    # Whiten mixed signals
    Xw = zca_whitening_matrix(Xc)

    # W
    W, obj_f, i = fastICA(Xw,  M, alpha=1)
    return W, obj_f, i
