import pandas as pd
from Dataset import Dataset

# Loading dataset
data_slice = snakemake.config['ICA_datasets'][snakemake.wildcards.dataset]['variables']
dataset = Dataset(snakemake.input.counts, data_slice)

# Loading components
components = pd.read_csv(
    snakemake.input.components,
    sep='\t', index_col=[0, 1, 2, 3]
)

# Removing unnecessary genes
keep_genes = components.index.tolist()
dataset.data = dataset.data.loc[keep_genes]

# Calculating dot product
dataset_projection = dataset.data.T.dot(components)

# Output results
dataset_projection.to_csv(snakemake.output.projection, sep='\t')
