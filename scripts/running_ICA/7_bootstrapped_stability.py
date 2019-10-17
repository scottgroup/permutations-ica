import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import seaborn as sns

"""
    Objective : attest the stability of the components in relation with the
    datasets

    Is a said component present in all iteration?
        Keeping only component that have a correlation > min_correlation
        Quantify the number of iterations where a component was found

    Is a said component stable?
        Quantify correlation of accepted bootstrapped components

"""

# Loading components from the model
components = pd.read_csv(
    snakemake.input.components, sep='\t', index_col=[0,1,2,3]
)

# Loading strapped components
boot_components = pd.read_csv(
    snakemake.input.boot_components, sep='\t', index_col=[0,1,2,3]
)

# Components - Boot-components association
comp2boot = dict()
boot_stability = dict()
for c in components.columns:
    comp2boot[c] = list()
    boot_stability[c] = list()


for c_boot in boot_components.columns:
    score = dict()
    for c in components.columns:
        r2, pval = pearsonr(components[c], boot_components[c_boot])
        score[c] = r2
    max_comp = max(score, key=lambda k: abs(score[k]))
    if abs(score[max_comp]) > snakemake.params.min_correlation:
        comp2boot[max_comp].append(c_boot)
        boot_stability[max_comp].append(score[max_comp])


data = [np.abs(boot_stability[c]) for c in boot_stability]
sns.violinplot(data=data, scale='count', inner='stick')
plt.savefig(snakemake.output.plot)
