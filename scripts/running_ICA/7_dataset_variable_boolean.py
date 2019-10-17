import pandas as pd
from Dataset import Dataset

# Loading dataset
data_slice = snakemake.config['ICA_datasets'][snakemake.wildcards.dataset]
dataset = Dataset(snakemake.input.counts, data_slice)

# Keeping index
dataset.data = dataset.data.T
dataset.data = dataset.data.drop(
    dataset.data.columns[:], axis=1
)
df = pd.DataFrame(index=dataset.data.index)

# Adding columns with first degree boolean
for levels, codes, name in zip(df.index.levels, df.index.codes, df.index.names):
    for level, i in zip(levels, range(len(levels))):
        df[level] = [1 if c == i else 0 for c in codes]

# Writing df to file
df.to_csv(snakemake.output.var_bool, sep='\t')
