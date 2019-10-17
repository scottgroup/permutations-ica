import pandas as pd
from scipy.stats import pearsonr

# Loading var_bool
var_bool = pd.read_csv(
    snakemake.input.var_bool, sep='\t', index_col=[0, 1, 2, 3, 4, 5]
)

# Loading projections along components
proj_df = pd.read_csv(
    snakemake.input.projection, sep='\t', index_col=[0, 1, 2, 3, 4, 5]
)

# Calculation pearsonr
df = pd.DataFrame(index=proj_df.columns, columns=var_bool.columns)
# For every component
for comp in proj_df.columns:

    # Calculated against all variables
    for var in var_bool.columns:
        r2, pval = pearsonr(proj_df[comp], var_bool[var])
        df.loc[comp, var] = pval

df.to_csv(snakemake.output.proj_corr, sep='\t')
