import json
import pandas as pd

# With parameters
sigma = int(snakemake.wildcards.sigma)

# Loading components
mean_df = pd.read_csv(
    snakemake.input.components_mean,
    sep='\t', index_col=[0, 1, 2, 3]
)

# Filtering genes
for col in mean_df.columns:
    grand_mean = mean_df[col].mean(axis=0)
    grand_std = mean_df[col].std(axis=0)

    filt_up = mean_df[col] < (grand_mean + sigma * grand_std)
    filt_down = mean_df[col] > (grand_mean - sigma * grand_std)
    filt = filt_up & filt_down

    mean_df[col].loc[filt] = 0

mean_df = mean_df[(mean_df.sum(axis=1) != 0)]

mean_df.to_csv(snakemake.output.filt_genes, sep='\t')
