import pandas as pd

hgnc_df = pd.read_csv(snakemake.input.hgnc, sep='\t', dtype='str')
hgnc_df.fillna('NaN', inplace=True)

gene_ids = [
    ('NCBI Gene ID', 'NCBI Gene ID(supplied by NCBI)'),
    ('Ensembl gene ID', 'Ensembl ID(supplied by Ensembl)')
]

df_list = list()
for idx, row in hgnc_df.iterrows():
    _line = [row['HGNC ID'], row['Approved symbol']]

    # Adding the annotation
    for g1, g2 in gene_ids:
        gene1 = row[g1]
        gene2 = row[g2]
        if gene1 == gene2 and gene1 != 'NaN':
            _line.append(gene1)
        elif gene1 == 'NaN' and gene2 != 'NaN':
            _line.append(gene2)
        elif gene1 != 'NaN':
            _line.append(gene1)
        else:
            _line.append('')

    df_list.append(_line)

df = pd.DataFrame(df_list, columns=['HGNC_ID', 'symbol', 'ncbi_id', 'ensembl_id'])
df.to_csv(snakemake.output.hgnc, sep='\t', index=None)
