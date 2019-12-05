import pandas as pd
import re
import snakemake.shell as shell

columns = [
    'seq', 'source', 'feature', 'start', 'end',
    'score', 'strand', 'frame', 'attributes'
]

# Loading data
gtf_file = snakemake.input.gtf
gtf_gene = shell(
    "cat {gtf_file} | awk '{{if ($3==\"gene\") print $0}}' ",
    iterable=True
)

# Format GTF file
gtf_lines = [line.split('\t') for line in gtf_gene if line[0] != '#']
df_gene = pd.DataFrame(gtf_lines, columns=columns)

df_gene['gene_id'] = df_gene['attributes'].str.extract("gene_id \"(.*?)\"")
df_gene['gene_biotype'] = df_gene['attributes'].str.extract("gene_biotype \"(.*?)\"")

df_gene = df_gene[['gene_id', 'gene_biotype']]
df_gene.to_csv(snakemake.output.biotype, sep='\t', index=False)
