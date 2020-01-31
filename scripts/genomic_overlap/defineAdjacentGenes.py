import pandas as pd
from snakemake import shell

# Loading parameters
annotation = snakemake.wildcards.annotation
gene = snakemake.wildcards.gene
bed = snakemake.input.bed
gtf = snakemake.input.gtf

# Intersecting bed with GTF
hits = shell(
    "bedtools intersect -wb -a {bed} -b {gtf}",
    iterable=True
)
df = pd.DataFrame(
    [hit.split('\t') for hit in hits]
)

# Extracting overlapping genes
re_st = "{key} \"(.*?)\""
if 'ENSG' in gene:
    gene_id = re_st.format(key="gene_id")
else:
    gene_id = re_st.format(key="db_xref")
df['gene'] = df[df.columns[-1]].str.extract(gene_id)
genes = df['gene'].unique()
if not 'ENSG' in gene:
    genes = [
        gene for gene in genes
        if isinstance(gene, str) and 'GeneID' in gene
    ]

with open(snakemake.output.tsv, 'w') as f:
    for gene in genes:
        f.write('\t'.join([gene, annotation]) + '\n')
