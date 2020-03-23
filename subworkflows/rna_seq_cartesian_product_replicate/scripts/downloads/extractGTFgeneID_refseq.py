import pandas as pd

columns = [
    'seqid', 'source', 'type', 'start',
    'end', 'score', 'strand', 'phase', 'attributes'
]

# Loading gtf into pandas
gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=5,
    index_col=False, header=None, names=columns
)

# Keeping only genes
gtf = gtf[gtf['type'] == 'gene']

# Dropping all columns except attributes
gtf.drop(
    [c for c in columns if c is not 'attributes'],
    axis=1, inplace=True
)

# Adding columns
gtf['symbol'] = gtf['attributes'].str.extract("gene_id \"(.*?)\"")
gtf['ncbi_id'] = gtf['attributes'].str.extract("db_xref \"GeneID:(.*?)\"")
gtf['HGNC_ID'] = gtf['attributes'].str.extract("\"HGNC:(HGNC:.*?)\"")

# Dropping attributes
gtf.drop('attributes', axis=1, inplace=True)

# Removing all genes having _1 in the gene id. Duplicates from X/Y pseudoautosomal regions
gtf = gtf[~gtf['symbol'].str.contains('_[0-9]+', regex=True)]

# Writing to CSV
gtf.to_csv(snakemake.output.refseq_gene, sep='\t', index=None, header=None)
