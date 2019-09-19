import json
import pandas as pd

# Loading components
components = pd.read_csv(
    snakemake.input.components,
    sep='\t', index_col=[0, 1, 2, 3]
)

# Loading correlated_components
with open(snakemake.input.corr_components, 'r') as f:
    corr_components = json.load(f)

# Creating a mean and std DataFrame for the components
mean_df = pd.DataFrame(index=components.index)
std_df = pd.DataFrame(index=components.index)
for comp, val in corr_components.items():
    mean_df[comp] = components[val].mean(axis=1)
    std_df[comp] = components[val].std(axis=1)

mean_df.to_csv(snakemake.output.components_mean, sep='\t')
std_df.to_csv(snakemake.output.components_std, sep='\t')
