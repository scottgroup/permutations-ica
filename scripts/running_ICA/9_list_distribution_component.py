import pandas as pd

# Loading data
proj_corr = pd.read_csv(
    snakemake.input.proj_corr, sep='\t', index_col=[0]
)
var_bool = pd.read_csv(
    snakemake.input.var_bool, sep='\t', usecols=[0, 1, 2, 3, 4, 5]
)

# Dictionary with variable relation
var_dict = dict()
for col in var_bool.columns:
    for _val in set(var_bool[col]):
        var_dict[_val] = col

# Dictionary with minimum and component
min_dict = proj_corr.idxmin(axis=1).to_dict()

# Writing to file
with open(snakemake.output.list_comp, 'w') as f:
    min_dict = proj_corr.idxmin(axis=1).to_dict()
    for k, val in min_dict.items():
        _str = str(k) +'\t' + var_dict[val] + '\n'
        f.write(_str)
