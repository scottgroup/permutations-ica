import pandas as pd

columns = [
    'seqid', 'source', 'type', 'start',
    'end', 'score', 'strand', 'phase', 'attributes'
]

# Loading gtf into pandas
gtf = pd.read_csv(
    snakemake.input.gtf,
    compression='gzip',
    sep='\t', skiprows=5,
    index_col=False, header=None, names=columns
)

# Keeping only transcript
gtf = gtf[gtf['type'] == 'transcript']

# Dropping all columns except attributes
gtf.drop(
    [c for c in columns if c is not 'attributes'],
    axis=1, inplace=True
)

# Adding columns
gtf['gene_id'] = gtf['attributes'].str.extract("gene_id \"(.*?)\"")
gtf['transcript_id'] = gtf['attributes'].str.extract("transcript_id \"(.*?)\"")

# Dropping attributes
gtf = gtf[['gene_id', 'transcript_id']]

# Writing to CSV
gtf.to_csv(snakemake.output.gene_map, sep='\t', index=None)
