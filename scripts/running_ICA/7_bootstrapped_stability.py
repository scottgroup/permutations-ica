import numpy as np
import pandas as pd
from scipy.stats import pearsonr


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
for c in components.columns:
    comp2boot[c] = list()



for c_boot in boot_components.columns:
    print('For component boot ', c_boot)
    score = dict()
    for c in components.columns:
        r2, pval = pearsonr(components[c], boot_components[c_boot])
        print('\t', r2)
        score[c] = r2
    max_comp = max(score, key=lambda k: abs(score[k]))
    comp2boot[max_comp].append(c_boot)

for c in comp2boot:
    print(comp2boot[c])
