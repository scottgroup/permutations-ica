import pandas as pd


def slice_dataset(df, data_slice):
    """ For the variables in data_slice, only keep specified values """
    df = df.transpose()
    slice_list = list()
    for idx in df.index.names:
        if idx in data_slice:
            slice_list.append(data_slice[idx])
        else:
            slice_list.append(slice(None))

    return df.loc[tuple(slice_list),].transpose()


# Loading the dataset
data = pd.read_csv(
    snakemake.input.dataset, sep='\t',
    index_col=[0, 1, 2, 3], header=[0, 1, 2, 3, 4, 5],
    na_values=[''],
    dtype={
        'HGNC_ID': 'object',
        'symbol': 'object',
        'ncbi_id': 'object',
        'ensembl_id': 'object',
    }
)

# Dropping non-full rows
data.dropna(axis=0, inplace=True)

# Removing unused variables
data_slice = snakemake.config['ICA_models'][snakemake.wildcards.ICAmodel]['variables']
data = slice_dataset(data, data_slice)

# Compression the multi-columns
data.columns = data.columns.map('|'.join).str.strip('|')

# Keep only one index
data.index = data.index.droplevel(['HGNC_ID', 'ncbi_id', 'ensembl_id'])
data.reset_index(inplace=True)

# Generating the sample table
samples = pd.DataFrame(
    [[c, *c.split('|')] for c in data.columns[1:]],
    columns=[
        'sample', 'quantifier', 'tissue',
        'dataset', 'trimmer', 'aligner', 'annotation'
    ]
)

# Data to file
data.to_csv(snakemake.output.counts, sep='\t', index=False)
samples.to_csv(snakemake.output.samples, sep='\t', index=False)
