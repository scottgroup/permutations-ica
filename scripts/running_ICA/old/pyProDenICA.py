import numpy as np
from rpy2.robjects.packages import STAP
import rpy2.robjects as robj
import os

with open('scripts/running_ICA/GFunc.r', 'r') as f:
    string = f.read()
GFunc = STAP(string, "GFunc")


with open('scripts/running_ICA/ProDenICA.r', 'r') as f:
    string = f.read()
ProDenICAfun = STAP(string, "ProDenICA")


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


def covarianceGit(x):
    mean = np.mean(x, axis=1, keepdims=True)
    n = np.shape(x)[1] - 1
    m = x - mean

    return (m.dot(m.T))/n

def whitenGit(x):
    # Calculate the covariance matrix
    coVarM = covarianceGit(x)

    # Single value decoposition
    U, S, V = np.linalg.svd(coVarM)

    # Calculate diagonal matrix of eigenvalues
    d = np.diag(1.0 / np.sqrt(S))

    # Calculate whitening matrix
    whiteM = np.dot(U, np.dot(d, U.T))

    # Project onto whitening matrix
    Xw = np.dot(whiteM, x)

    return Xw, whiteM


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


def ProDenICA(X,
              k=None,
              W0=None,
              whiten=False,
              maxit=500,
              thresh=1e-7,
              restarts=0,
              trace=False,
              Gfunc=ProDenICAfun.GPois,
              eps_rank=1e-7):
    """ Running the R implementation """
    print('This is X')
    print(X)
    print(type(X))

    data = numpy2r(X)

    ProDenICAfun.ProDenICA(
        x=data, trace=True, maxit=10, Gfunc=Gfunc
        # k=k, W0=W0, whiten=whiten,,
        # maxit=maxit, thresh=thresh, restarts=restarts,
        # trace=trace# , # Gfunc=Gfunc, eps_rank=eps_rank
    )


def pyProDenICA(X,
                k=None,
                W0=None,
                whiten=False,
                maxit=500,
                thresh=1e-7,
                restarts=0,
                trace=False,
                Gfunc=GFunc.G1,
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
    # X -= np.mean(X, axis=0)

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
        # X, _ = whitenGit(X)

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

    gpS_file = 'gpS_file.txt'
    gS_file = 'gS_file.txt'

    for n_it in range(maxit):
        # Running something similar as previous
        s = np.dot(X, W0)
        for i in range(len(flist)):
            flist[i] = parse_GFunc(Gfunc(numpy2r(s[:, i])))

        crit0 = np.mean([d['Gs'] for d in flist0])
        gS = np.array([d['gs'] for d in flist0])
        gpS = np.array([d['gps'] for d in flist0])

        with open(gS_file, 'a') as f:
            f.write(str([list(x) for x in gS]) + '\n')

        with open(gpS_file, 'a') as f:
            f.write(str([list(x) for x in gpS]) + '\n')
        # print(crit0)
        # print(gS)
        # print(gpS)

        t1 = np.dot(X.T, gS.T/n)
        t2 = np.mean(gpS, axis=1)
        W1 = t1 - W0/(1/t2)
        W1 = ICAorthW(W1)
        nw = amari(W0, W1)

        if trace:
            print('Iteration {n_it} G {crit0} crit {nw}'.format(n_it=n_it, crit0=crit0, nw=nw))

        W0 = W1

        # Reach stability
        if nw < thresh:
            return W0, s, crit0, nw, n_it

    return 0, 0, crit0, nw, n_it
