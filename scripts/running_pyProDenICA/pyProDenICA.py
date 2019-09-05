import numpy as np
from rpy2.robjects.packages import STAP
import rpy2.robjects as robj
import os

print(os.getcwd())


with open('scripts/running_pyProDenICA/GFunc.r', 'r') as f:
    string = f.read()
GFunc = STAP(string, "GFunc")


def parse_GFunc(gfunc):
    return dict(zip(gfunc.names, map(list, list(gfunc))))


def ylim_scale(ylim, scale=0):
    scale2 = np.diff(ylim)

    if scale2 < scale:
        return np.repeat(np.mean(ylim), 2) + ((ylim - np.mean(ylim)) * scale)/scale2
    else:
        return ylim


def numpy2r(mat):
    if len(mat.shape) == 2:
        nr, nc = mat.shape
        xvec = robj.FloatVector(mat.transpose().reshape((mat.size)))
        return robj.r.matrix(xvec, nrow=nr, ncol=nc)
    else:
        return robj.FloatVector(mat)


def ICAorthW(W):
    u, v, sh = np.linalg.svd(W, full_matrices=False)
    return np.dot(u, sh)


def amari(V, W, orth=True):
    """ """
    if orth:
        A = abs(np.dot(V.T, W))
    else:
        A = abs(np.linalg.solve(V, W))

    rmax = np.max(A, axis=1)
    cmax = np.max(A, axis=0)
    rsum = np.sum(A, axis=1)
    csum = np.sum(A, axis=0)

    return (sum(rsum/rmax - 1) + sum(csum/cmax - 1))/(2 * np.size(A, 0))


def zca_whitening_matrix(X):
    """
    Function to compute ZCA whitening matrix (aka Mahalanobis whitening).
    INPUT:  X: [M x N] matrix.
        Rows: Variables
        Columns: Observations
    OUTPUT: ZCAMatrix: [M x M] matrix
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
    ZCAMatrix = np.dot(U, np.dot(np.diag(1.0/np.sqrt(S + epsilon)), U.T))  # [M x M]
    return ZCAMatrix


def pyProDenICA(X,
                k=None,
                W0=None,
                whiten=False,
                maxit=500,
                thresh=1e-7,
                restarts=0,
                trace=False,
                Gfunc=GFunc.GPois,
                eps_rank=1e-7):
    """Performs ProDenICA analysis

    Parameters
    _______________

    X : (,) array
    """
    # Data dimensions
    n, p = np.shape(X)
    if k is None:
        k = p

    # Centering X
    X -= np.mean(X, axis=0)

    # Whitening data
    if whiten:
        # SVD
        print(X.shape)
        # u, v, sh = np.linalg.svd(X)
        # condnum = v/v[0]
        # good = condnum > eps_rank
        # rank = np.sum(good)
        #
        # # Reducing k to match rank
        # if k > rank:
        #     print("rank of x is {rank}; k reduced from {k} to {rank}".format(
        #         k=k, rank=rank
        #     ))
        #     k = rank
        # X = np.sqrt(n) * v[:, good]

        # ZCA whitening
        zca_whitener = zca_whitening_matrix(X)
        X = np.dot(zca_whitener, X)

        print('After whitening ', X.shape)
    else:
        whitener = None

    # If no W0, generate starting point
    if W0 is None:
        W0 = np.random.normal(p, k, (p, k))
    else:
        k = np.shape(W0)[1]
    W0 = ICAorthW(W0)

    # Initialisation
    Gs = np.zeros((n, k))
    gS = Gs
    gpS = Gs

    s = np.dot(X, W0)

    flist = list(np.arange(0, k))
    for i in range(len(flist)):
        flist[i] = parse_GFunc(Gfunc(numpy2r(s[:, i])))

    flist0 = flist
    crit0 = np.mean([d['Gs'] for d in flist0])

    # ! Restarts not yet implemented
    for i in range(restarts):
        W1 = np.random.normal(p, k, (p, k))
        W1 = ICAorthW(W1)
        s = np.dot(X, W1)
        for j in range(k):
            flist[j] = parse_GFunc(Gfunc(numpy2r(s[:, j])))
        crit = np.mean([d['Gs'] for d in flist0])

        if trace:
            print('old crit {crit0} new crit {crit}'.format(crit0=crit0, crit=crit))

        if crit > crit0:
            crit0 = crit
            W0 = W1
            flist0 = flist

    # Main loop
    nw = 10

    for n_it in range(maxit):
        # Running something similar as previous
        s = np.dot(X, W0)
        for i in range(len(flist)):
            flist[i] = parse_GFunc(Gfunc(numpy2r(s[:, i])))

        crit0 = np.mean([d['Gs'] for d in flist0])
        gS = np.array([d['gs'] for d in flist0])
        gpS = np.array([d['gps'] for d in flist0])

        t1 = np.dot(X.T, gS.T/n)  # No tranpose? Due to gS being linear?
        t2 = np.mean(gpS, axis=1)
        W1 = t1 - W0/(1/t2)
        W1 = ICAorthW(W1)
        nw = amari(W0, W1)

        if trace:
            print('Iteration {n_it} G {crit0} crit {nw}'.format(n_it=n_it, crit0=crit0, nw=nw))

        W0 = W1

        # Reach stability
        if nw < thresh:
            return W0, crit0, nw, n_it

    return 0, crit0, nw, n_it
