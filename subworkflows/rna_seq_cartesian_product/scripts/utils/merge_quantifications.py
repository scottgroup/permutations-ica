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

# # Updating HGNC_df to add index for genes without HGNC ID
# for annotation in snakemake.input.quants:
#
#     if 'ensembl' in annotation:
#         index_col = [0, 1]
#     elif 'refseq' in annotation:
#         index_col = [0, 1, 2]
#
#     quant_df = pd.read_csv(
#         annotation, sep='\t',
#         index_col=index_col, header=[0, 1, 2, 3, 4, 5], dtype='str'
#     )
#     quant_df.index = quant_df.index.astype('object')
#
#     if 'ensembl' in annotation:
#         known_ENGS = hgnc_df['ensembl_id'].unique().tolist()
#         to_add = quant_df[
#             ~quant_df.index.get_level_values('gene').isin(known_ENGS)
#         ]
#         to_add = to_add.reset_index()[['gene', 'symbol']]
#
#         # Formatting DataFrame
#         to_add = pd.DataFrame(
#             [("", symbol, "", ens_id) for ens_id, symbol in to_add.values],
#             columns=["HGNC_ID", "symbol", "ncbi_id", "ensembl_id"]
#         )
#     elif 'refseq' in annotation:
#         known_NCBI = [int(c) for c in hgnc_df['ncbi_id'].unique().tolist() if c is not np.nan and c is not '']
#         to_add = quant_df[
#             ~quant_df.index.get_level_values('gene').isin(known_NCBI)
#         ]
#         to_add = to_add.reset_index()[['gene', 'symbol']]
#
#         # Formatting DataFrame
#         to_add = pd.DataFrame(
#             [("", symbol, ncbi_id, "") for ncbi_id, symbol in to_add.values],
#             columns=["HGNC_ID", "symbol", "ncbi_id", "ensembl_id"]
#         )
#         to_add.index = to_add.index.astype('object')
#     hgnc_df = pd.concat([hgnc_df, to_add])


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
hgnc_df.to_csv(snakemake.output.out_file, sep='\t')

hgnc_df.dropna(axis=0, inplace=True)
