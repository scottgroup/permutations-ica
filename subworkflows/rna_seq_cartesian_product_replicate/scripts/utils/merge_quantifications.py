import numpy as np
import pandas as pd


def get_index_col(annotation):
    if 'ensembl' in annotation:
        return [0, 1]
    elif 'refseq' in annotation:
        return [0, 1, 2]

def get_merge_field(annotation):
    if 'ensembl' in annotation:
        return 'gene', 'ensembl_id'
    elif 'refseq' in annotation:
        return 'HGNC_ID', 'HGNC_ID'


hgnc_df = pd.read_csv(snakemake.input.hgnc, sep='\t', dtype='str')

# Adding data
for annotation in snakemake.input.quants:
    right_on, left_on = get_merge_field(annotation)
    index_col = get_index_col(annotation)

    quant_df = pd.read_csv(
        annotation, sep='\t',
        index_col=index_col, header=[0, 1, 2, 3, 4, 5], dtype='str'
    )

    hgnc_df = hgnc_df.merge(
        right=quant_df, how='left',
        right_on=right_on, left_on=left_on
    )

# Keeping only HGNC
hgnc_df = hgnc_df[hgnc_df['HGNC_ID'].str.contains('HGNC')]

hgnc_df.set_index(['HGNC_ID', 'symbol', 'ncbi_id', 'ensembl_id'], inplace=True)

hgnc_df.columns = pd.MultiIndex.from_tuples(
    hgnc_df.columns,
    names=['quantifier', 'tissue', 'data_id', 'trimmer', 'aligner', 'annotation']
)

if snakemake.wildcards.isNaN == 'noNaN':
    hgnc_df.fillna(0, axis=1, inplace=True)
else:
    hgnc_df.dropna(axis=0, inplace=True)

hgnc_df.to_csv(snakemake.output.out_file, sep='\t')
