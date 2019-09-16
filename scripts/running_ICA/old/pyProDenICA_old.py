import numpy as np


def G1(s, a=1):
    return {
        "Gs": np.log(np.cosh(a * s))/a,
        "gs": np.tanh(a * s),
        "gps": a * (1 - np.tanh(a * s)**2)
    }


def ylim_scale(ylim, scale=0):
    scale2 = np.diff(ylim)

    if scale2 < scale:
        return np.repeat(np.mean(ylim), 2) + ((ylim - np.mean(ylim)) * scale)/scale2
    else:
        return ylim


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
                Gfunc=G1,
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
        # ZCA whitening
        zca_whitener = zca_whitening_matrix(X)
        X = np.dot(zca_whitener, X)
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
        flist[i] = G1(s[:, i])
        print(flist[i])

    flist0 = flist
    crit0 = np.mean([d['Gs'] for d in flist0])

    # Restarts
    for i in range(restarts):
        W1 = np.random.normal(p, k, (p, k))
        W1 = ICAorthW(W1)
        s = np.dot(X, W1)
        for j in range(k):
            flist[j] = G1(s[:, j])
        crit = np.mean([d['Gs'] for d in flist0])

        if trace:
            print('old crit {crit0} new crit {crit}'.format(crit0=crit0, crit=crit))

        if crit > crit0:
            crit0 = crit
            W0 = W1
            flist0 = flist

    # Main loop
    for n_it in range(maxit):
        # Running something similar as previous
        print(X)
        print('X shape')
        print(X.shape)
        s = np.dot(X, W0)
        for i in range(len(flist)):
            flist[i] = G1(s[:, i])

        crit0 = np.mean([d['Gs'] for d in flist0])
        gS = np.array([d['gs'] for d in flist0])
        gpS = np.array([d['gps'] for d in flist0])

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
            print('W0 shape')
            print(W0.shape)
            print('s shape')
            print(s.shape)
            return W0, crit0, nw, n_it

    return 0, crit0, nw, n_it
