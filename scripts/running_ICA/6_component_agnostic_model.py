import json
import pandas as pd
from scipy.stats import pearsonr


def find_correlation(components, col, corr_dict):
    """ Creates dictionnary with correlated components """
    top_hit, top_key = 0, ''

    # For every already defined clusters
    for key in corr_dict.keys():

        # Calculate correlation with every element of the cluster
        pear = list()
        for corr_col in corr_dict[key]:
            try:
                r2, pval = pearsonr(components[col], components[corr_col])
            except ValueError:
                print('col ', col)
                print('corr_col ', corr_col)
                print(components[col])
                print(components[corr_col])
                raise
            pear.append(r2)
        sum_pear = sum([abs(x) for x in pear])/len(pear)
        if sum_pear > top_hit:
            top_hit = sum_pear
            top_key = key

    return top_hit, top_key


def flip_components(components, corr_dict):
    """ Returns a list of components that have opposite direction """
    to_flip = list()

    for key in corr_dict.keys():
        base_comp = corr_dict[key][0]
        for comp in corr_dict[key][1:]:
            r2, pval = pearsonr(components[base_comp], components[comp])
            if r2 < 0:
                to_flip.append(comp)
    return to_flip


# Loading wildcards
M_min = int(snakemake.wildcards.min)
M_max = int(snakemake.wildcards.max)
n = int(snakemake.wildcards.n)
components_file = snakemake.input

# Loading all components into a dataframe
components = pd.DataFrame()
for file, M in zip(components_file, range(M_min, M_max+1)):
    _component = pd.read_csv(
        file,
        sep='\t', index_col=[0, 1, 2, 3]
    )
    _component.columns = [col + '_M' + str(M) for col in _component.columns]
    components = pd.concat([components, _component], axis=1, sort=False)

# Sorting the components into clusters, based on a threshold
threshold = snakemake.params.threshold
n = int(snakemake.wildcards.n)

# Solving though iterations
corr_co = dict()
for col in components.columns.tolist():

    # Create first dictionary entry
    if len(corr_co) == 0:
        corr_co[len(corr_co)] = [col]
    else:
        hit, key = find_correlation(components, col, corr_co)
        if hit > threshold:
            corr_co[key].append(col)
        else:
            corr_co[len(corr_co)] = [col]

# Flipping components
to_flip = flip_components(components, corr_co)
for col in to_flip:
    components[col] *= -1


# Output of components
components.to_csv(snakemake.output.components, sep='\t')

# Output of correlated components
with open(snakemake.output.corr_components, 'w') as f:
    json.dump(corr_co, f)
