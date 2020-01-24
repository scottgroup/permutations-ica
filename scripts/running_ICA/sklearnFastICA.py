import numpy as np
from sklearn.decomposition import FastICA
from sklearn.exceptions import ConvergenceWarning


def sklearnFastICA(data, M, max_it, tolerance):
    """
    Using the FastICA implementation from sklearn
    """
    # Constructing the ICA object
    ICA = FastICA(
        n_components=M, whiten=True,
        algorithm='parallel', fun='cube',
        max_iter=max_it, tol=tolerance,
    )

    try:
        ICA.fit(np.array(data.T))

        if ICA.n_iter_ == max_it:
            components = ICA.components_.T
            return components
        else:
            return None

    except ValueError as err:
        print("Value error: {err}".format(err=err))
        return None
