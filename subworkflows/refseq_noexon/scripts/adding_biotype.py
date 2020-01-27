import pandas as pd
from snakemake import shell

def loading_csv(path):
    return pd.read_csv(
        path, sep='\t', skiprows=[0,1,2,3,4],
        names=[
            'seqname', 'source', 'feature', 'start',
            'end', 'score', 'strand', 'frame', 'attributes'
        ]
    )

# RefSeq biotypes
df = loading_csv(snakemake.input.gtf_refseq)
df = df[df['feature'] == 'gene']
df['GeneID'] = df['attributes'].str.extract("GeneID:(.*?)\"")
df['biotype'] = df['attributes'].str.extract("gene_biotype \"(.*?)\"")
df = df[['GeneID', 'biotype']]
refseq_biotype_dict = df.set_index('GeneID')['biotype'].to_dict()

# Ensembl biotypes
df = loading_csv(snakemake.input.gtf_ens98)
df = df[df['feature'] == 'gene']
df['GeneID'] = df['attributes'].str.extract("gene_id \"(.*?)\"")
df['biotype'] = df['attributes'].str.extract("gene_biotype \"(.*?)\"")
df = df[['GeneID', 'biotype']]
ensembl_biotype_dict = df.set_index('GeneID')['biotype'].to_dict()

# Loading the gene counts
data = pd.read_csv(
    snakemake.input.hgnc_noRNA, sep='\t',
    dtype={
        'HGNC_ID': 'object',
        'symbol': 'object',
        'ncbi_id': 'object',
        'ensembl_id': 'object',
    }
)
data['refseq_biotype'] = data['ncbi_id'].map(refseq_biotype_dict)
data['ensembl_biotype'] = data['ensembl_id'].map(ensembl_biotype_dict)
data.to_csv(snakemake.output.result, sep='\t', index=False)
