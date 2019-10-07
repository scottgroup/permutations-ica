import numpy as np
from rpy2.robjects.packages import STAP
import rpy2.robjects as robj

# Loading R ConsICA script
with open('scripts/running_ICA/ICAmethods/consICA.R', 'r') as f:
    string = f.read()
runICA = STAP(string, "runICA")


def numpy2r(mat):
    if len(mat.shape) == 2:
        nr, nc = mat.shape
        xvec = robj.FloatVector(mat.transpose().reshape((mat.size)))
        return robj.r.matrix(xvec, nrow=nr, ncol=nc)
    else:
        return robj.FloatVector(mat)


# running ConsICA
def consICA(X, ncomp, ntry, ncores):
    """ """
    data = numpy2r(np.array(X))
    runICA.runICA(X=data, ncomp=ncomp, ntry=ntry, ncores=ncores)
