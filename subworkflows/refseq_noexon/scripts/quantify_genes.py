import pandas as pd

# Reading HGNC
hgnc = pd.read_csv(
    snakemake.input.hgnc, sep='\t',
    dtype={
        'HGNC_ID': 'object',
        'symbol': 'object',
        'ncbi_id': 'object',
        'ensembl_id': 'object',
    }
).dropna()

# Reading gene counts
data = pd.read_csv(snakemake.input.GeneID_numbers, sep='\t')

# Keeping genes with unique count
data = data[data['N'] == 1]
genes = [str(gene) for gene in data['GeneID'].tolist()]

# Slicing HGNC list
hgnc_noRNA = hgnc[hgnc['ncbi_id'].isin(genes)]
hgnc_noRNA.to_csv(snakemake.output.hgnc_noRNA, sep='\t', index=None)
