import numpy as np
from sklearn.decomposition import FastICA
from sklearn.exceptions import ConvergenceWarning

# import warnings
# warnings.filterwarnings("error")

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

    # try:
    #     # Processing the dataset
    #     ICA.fit(np.array(data.T)) # n_samples, n_features
    #     # print(dir(ICA))
    #     print(ICA.n_iter_)
    #
    #     # Extracting the components
    #     components = ICA.components_.T # n_features, n_components
    #     return components
    #
    # except ConvergenceWarning:
    #     return None

    ICA.fit(np.array(data.T))

    print(ICA.n_iter_)

    if ICA.n_iter_ == max_it:
        components = ICA.components_.T
        return components
    else:
        return None
