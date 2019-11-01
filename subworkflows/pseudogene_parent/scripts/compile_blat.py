import pandas as pd

cols = ['query', 'subject', 'e-value', 'bit_score']
df = pd.DataFrame(columns=cols)

for blat in snakemake.input.blats:
    data = pd.read_csv(blat, sep='\t', names=[
        'query', 'subject', 'identity', 'alignment', 'mismatches', 'gap',
        'q_start', 'q_end', 's_start', 's_end', 'e-value', 'bit_score'
    ])

    _df = data.sort_values('e-value', ascending=True).drop_duplicates('query')
    df = df.append(_df[cols])

# Adding the gene column
df_gene = pd.read_csv(snakemake.input.gene_map, sep='\t', index_col=1)
df['subjectGene'] = df['subject'].map(df_gene['gene_id'].to_dict())

df.to_csv(snakemake.output.score, sep='\t', index=None)
