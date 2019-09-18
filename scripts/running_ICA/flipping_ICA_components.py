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
            r2, pval = pearsonr(components[col], components[corr_col])
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


components = pd.read_csv(
    snakemake.input.components,
    sep='\t', index_col=[0, 1, 2, 3]
)

# Sorting the components into clusters, based on a threshold
threshold = snakemake.params.threshold
minimal_occurence = snakemake.params.minimal_occurence
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

# Dropping
to_drop = list()
for key in corr_co:
    if len(corr_co[key]) < n * minimal_occurence:
        to_drop.append(key)
for dropping in to_drop:
    del corr_co[dropping]

# Flipping components
to_flip = flip_components(components, corr_co)
for col in to_flip:
    components[col] *= -1

# Keeping components that have N elements, all from different iterations
good_corr_co = dict()
good_comps = list()
for key, cluster in corr_co.items():
    if len(cluster) == n:
        iterations = [it.split()[0][2:] for it in cluster]
        if len(iterations) == len(set(iterations)):
            good_corr_co[key] = cluster
            good_comps.extend(cluster)

components = components[good_comps]

# Output of correlated components
with open(snakemake.output.corr_components, 'w') as f:
    json.dump(good_corr_co, f)

# Output of flipped components
components.to_csv(snakemake.output.components, sep='\t')
