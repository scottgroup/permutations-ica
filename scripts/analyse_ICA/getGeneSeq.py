import pandas as pd
from snakemake import shell

# Get gene informations
gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=[0,1,2,3,4],
    names=[
        'seqname', 'source', 'feature', 'start',
        'end', 'score', 'strand', 'frame', 'attributes'
    ]
)
gtf['gene_id'] = gtf['attributes'].str.extract("gene_id \"(.*?)\"")
gtf = gtf[gtf['feature'] == 'gene']
gene_gtf = gtf[gtf['gene_id'] == snakemake.wildcards.gene]
for idx, row in gene_gtf.iterrows():
    chr, strand = row['seqname'], row['strand']
    start, end = row['start'], row['end']

# Running samtools
shell(
    "samtools faidx {snakemake.input.genome} {chr}:{start}-{end} > {snakemake.output.gene}"
)
