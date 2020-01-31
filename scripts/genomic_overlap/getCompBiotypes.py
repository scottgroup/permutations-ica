import pandas as pd

def reading_file(path):
    up, down = list(), list()
    with open(path, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line == '>Positive genes':
                reading_up = True
            elif line == '>Negative genes':
                reading_up = False
            else:
                if reading_up:
                    up.append(line)
                else:
                    down.append(line)
    return up, down


gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=[0,1,2,3,4], usecols=[8],
    names=[
        'seqname', 'source', 'feature', 'start',
        'end', 'score', 'strand', 'frame', 'attributes'
    ]
)

gtf['gene_id'] = gtf['attributes'].str.extract("gene_id \"(.*?)\"")
gtf['biotype'] = gtf['attributes'].str.extract("gene_biotype \"(.*?)\"")
gtf = gtf[['gene_id', 'biotype']]
gtf.drop_duplicates(inplace=True)

# Getting the gene lists
up, down = reading_file(snakemake.input.gene_list)
s1 = gtf[gtf['gene_id'].isin(up)].groupby('biotype').count()
s2 = gtf[gtf['gene_id'].isin(down)].groupby('biotype').count()
biotypes = set([k for s in [s1, s2] for k in s.index.tolist()])
df = pd.DataFrame(index=biotypes)
df['up'] = df.index.map(s1.to_dict()['gene_id'])
df['down'] = df.index.map(s2.to_dict()['gene_id'])
df.to_csv(snakemake.output.gene_biotypes, sep='\t')
