import numpy as np
from rpy2.robjects.packages import STAP
import rpy2.robjects as robj

# Loading R ConsICA script
with open('scripts/running_ICA/ICAmethods/consICA.R', 'r') as f:
    string = f.read()
runICA = STAP(string, "runICA")

# Loading dataset

# running ConsICA
def consICA(X, ncomp, ntry, show_every, ncores):

    runICA.runICA()
