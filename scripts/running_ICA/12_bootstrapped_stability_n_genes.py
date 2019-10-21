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

def keeping_sigma_genes(df, sigma):
    _df = df.copy(deep=True)
    for col in _df.columns:
        grand_mean = _df[col].mean(axis=0)
        grand_std = _df[col].std(axis=0)

        filt_up = _df[col] < (grand_mean + sigma * grand_std)
        filt_down = _df[col] > (grand_mean - sigma * grand_std)
        filt = filt_up & filt_down

        _df[col].loc[filt] = 0

    return _df[(_df.sum(axis=1) != 0)]


def extract_genes(c):
    return [
        c[c > 0].index.get_level_values('symbol').tolist(),
        c[c < 0].index.get_level_values('symbol').tolist()
    ]

def get_gene_overlap(upref, downref, uptest, downtest):
    overlap_up = [g for g in upref if g in uptest]
    overlap_down = [g for g in downref if g in downtest]

    # overflow_up = [g for g in uptest if g not in upref]
    # overflow_down = [g for g in downtest if g not in downref]

    return max([
        len(overlap_up / len(upref)),
        len(overlap_down / len(downref)),
    ])
    # return (len(overlap_up) + len(overlap_down)) / (len(upref) + len(downref))


def quant_similar_genes(cref, ctest):
    upref, downref = extract_genes(cref)
    uptest, downtest = extract_genes(ctest)

    quant = get_gene_overlap(upref, downref, uptest, downtest)
    quant_inv = get_gene_overlap(upref, downref, downtest, uptest)

    return max([quant, quant_inv])


# Loading components from the model
components = pd.read_csv(
    snakemake.input.components, sep='\t', index_col=[0,1,2,3]
)

# Loading strapped components
boot_components = pd.read_csv(
    snakemake.input.boot_components, sep='\t', index_col=[0,1,2,3]
)

# Loading sigma
sigma = snakemake.params.sigma

print(components)
components = keeping_sigma_genes(components, sigma)
print(components)

boot_components = keeping_sigma_genes(boot_components, sigma)
print(boot_components)

# Components - Boot-components association
comp2boot = dict()
boot_stability = dict()
for c in components.columns:
    comp2boot[c] = list()
    boot_stability[c] = list()


for c_boot in boot_components.columns:
    score = dict()
    for c in components.columns:
        _score = quant_similar_genes(components[c], boot_components[c_boot])
        score[c] = _score
    max_comp = max(score, key=lambda k: abs(score[k]))
    if abs(score[max_comp]) > snakemake.params.min_genes:
        comp2boot[max_comp].append(c_boot)
        boot_stability[max_comp].append(score[max_comp])

data = [np.abs(boot_stability[c]) for c in boot_stability]

sns.violinplot(data=data, scale='count', inner='stick')
plt.savefig(snakemake.output.plot)
